//! GTDB-Tk-style log lines (`[YYYY-MM-DD HH:MM:SS] INFO: message`) written
//! to stderr and to `{out_dir}/{prefix}.log`.

use chrono::Local;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;
use std::sync::Mutex;

pub(crate) struct ClassifyLog {
    log: Mutex<BufWriter<File>>,
    warnings: Mutex<BufWriter<File>>,
}

impl ClassifyLog {
    pub(crate) fn create(out_dir: &Path, prefix: &str) -> io::Result<Self> {
        Ok(Self {
            log: Mutex::new(BufWriter::new(File::create(
                out_dir.join(format!("{prefix}.log")),
            )?)),
            warnings: Mutex::new(BufWriter::new(File::create(
                out_dir.join(format!("{prefix}.warnings.log")),
            )?)),
        })
    }

    fn line(level: &str, message: &str) -> String {
        format!(
            "[{}] {level}: {message}",
            Local::now().format("%Y-%m-%d %H:%M:%S")
        )
    }

    pub(crate) fn info(&self, message: impl AsRef<str>) {
        let line = Self::line("INFO", message.as_ref());
        eprintln!("{line}");
        let mut log = self.log.lock().expect("log mutex poisoned");
        let _ = writeln!(log, "{line}");
        let _ = log.flush();
    }

    pub(crate) fn warn(&self, message: impl AsRef<str>) {
        let line = Self::line("WARNING", message.as_ref());
        eprintln!("{line}");
        for sink in [&self.log, &self.warnings] {
            let mut sink = sink.lock().expect("log mutex poisoned");
            let _ = writeln!(sink, "{line}");
            let _ = sink.flush();
        }
    }
}
