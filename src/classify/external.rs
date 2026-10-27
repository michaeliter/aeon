//! Locating and running the external programs used by aeon classify.

use std::env;
use std::error::Error;
use std::path::{Path, PathBuf};
use std::process::Command;

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

/// Run a command to completion, turning a non-zero exit into an error that
/// carries the command line and its stderr.
pub(crate) fn run(command: &mut Command) -> Result<(), Box<dyn Error>> {
    let output = command.output().map_err(|error| {
        format!(
            "failed to start {}: {error}",
            command.get_program().to_string_lossy()
        )
    })?;
    if output.status.success() {
        return Ok(());
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
