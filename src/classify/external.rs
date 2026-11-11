//! Locating and running the external programs used by aeon classify.

use std::env;
use std::error::Error;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output, Stdio};
use std::thread;

/// Resolve `program` (a name or path) to an executable, searching PATH.
pub(crate) fn find_program(program: &str) -> Result<PathBuf, Box<dyn Error>> {
    let candidate = Path::new(program);
    if candidate.components().count() > 1 {
        return if candidate.is_file() {
            Ok(candidate.to_path_buf())
        } else {
            Err(format!("{program} does not exist").into())
        };
    }
    env::var_os("PATH")
        .into_iter()
        .flat_map(|paths| env::split_paths(&paths).collect::<Vec<_>>())
        .map(|dir| dir.join(program))
        .find(|path| path.is_file())
        .ok_or_else(|| {
            format!(
                "{program} was not found on PATH; install it (see environment.yml) or pass its path"
            )
            .into()
        })
}

/// `program` from `bin_dir` when given, else from PATH.
pub(crate) fn find_tool(bin_dir: Option<&str>, program: &str) -> Result<PathBuf, Box<dyn Error>> {
    match bin_dir {
        Some(dir) => find_program(&Path::new(dir).join(program).to_string_lossy()),
        None => find_program(program),
    }
}

/// Run a command to completion, turning a non-zero exit into an error that
/// carries the command line and its stderr.
pub(crate) fn run(command: &mut Command) -> Result<(), Box<dyn Error>> {
    run_output(command).map(|_| ())
}

/// `run`, returning the command's captured output.
pub(crate) fn run_output(command: &mut Command) -> Result<Output, Box<dyn Error>> {
    let output = command
        .output()
        .map_err(|error| start_error(command, error))?;
    check(command, output)
}

/// `run` with `feed` writing the command's standard input.
pub(crate) fn run_with_input(
    command: &mut Command,
    feed: impl FnOnce(&mut dyn Write) -> Result<(), Box<dyn Error>>,
) -> Result<(), Box<dyn Error>> {
    let mut child = command
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|error| start_error(command, error))?;
    let mut stdin = BufWriter::new(child.stdin.take().expect("piped stdin"));
    // Feed from this thread while another collects the output, so neither
    // side blocks on a full pipe; closing stdin ends the program's input.
    let (fed, output) = thread::scope(|scope| {
        let output = scope.spawn(|| child.wait_with_output());
        let fed = feed(&mut stdin).and_then(|()| Ok(stdin.flush()?));
        drop(stdin);
        (fed, output.join().expect("output thread panicked"))
    });
    let output = output?;
    // A program that fails early closes its stdin, so report its exit first.
    let output = check(command, output)?;
    fed.map(|()| drop(output))
}

fn start_error(command: &Command, error: std::io::Error) -> Box<dyn Error> {
    format!(
        "failed to start {}: {error}",
        command.get_program().to_string_lossy()
    )
    .into()
}

fn check(command: &Command, output: Output) -> Result<Output, Box<dyn Error>> {
    if output.status.success() {
        return Ok(output);
    }
    let arguments = command
        .get_args()
        .map(|argument| argument.to_string_lossy())
        .collect::<Vec<_>>()
        .join(" ");
    Err(format!(
        "{} {arguments} exited with {}: {}",
        command.get_program().to_string_lossy(),
        output.status,
        String::from_utf8_lossy(&output.stderr).trim()
    )
    .into())
}
