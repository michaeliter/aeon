//! Index warm-up. Search memory-maps an index and follows graph edges in
//! random order, so on a cold page cache every step waits on a random read.
//! Reading the file front to back first is far cheaper: on GTDB-sized
//! indexes on HPC storage a cold search ran about 240 times slower than a
//! warm one, and the sequential read took about 3 s. Files are read in a
//! background thread so callers can overlap them with other work, and are
//! skipped when they would not fit in the memory available to this process.

use std::collections::HashMap;
use std::fmt;
use std::fs::{self, File};
use std::io::{self, Read};
use std::path::{Path, PathBuf};
use std::sync::Mutex;
use std::sync::mpsc::{self, Receiver};
use std::thread;
use std::time::Instant;

use crate::metadata::index_path;

/// What warming one index did.
#[derive(Debug, Clone, PartialEq)]
pub(crate) enum Preload {
    Loaded {
        path: String,
        bytes: u64,
        seconds: f64,
    },
    /// The file would not fit in the memory available to this process.
    Skipped {
        path: String,
        bytes: u64,
        available: u64,
    },
    /// Warm-up was turned off (`--no-preload`).
    Disabled,
}

impl fmt::Display for Preload {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Loaded {
                path,
                bytes,
                seconds,
            } => write!(
                f,
                "Preloaded {path} ({}) in {seconds:.1}s",
                gigabytes(*bytes)
            ),
            Self::Skipped {
                path,
                bytes,
                available,
            } => write!(
                f,
                "Not preloading {path} ({}): only {} of memory is available, so the search reads it on demand",
                gigabytes(*bytes),
                gigabytes(*available)
            ),
            Self::Disabled => write!(f, "Index preloading is off"),
        }
    }
}

fn gigabytes(bytes: u64) -> String {
    format!("{:.2} GB", bytes as f64 / 1e9)
}

enum Pending {
    Running(Receiver<io::Result<Preload>>),
    Done(Preload),
}

/// Background warm-ups of index files, keyed by index prefix.
pub(crate) struct Warmup {
    pending: Mutex<HashMap<String, Pending>>,
}

impl Warmup {
    /// Start reading the index files of `prefixes`, one after another in a
    /// background thread, so the first is ready soonest and the reads do not
    /// compete for I/O. Indexes are admitted in order while they fit, together,
    /// in the memory available to this process; the rest are skipped. With
    /// `enabled` false nothing is read.
    pub(crate) fn start(prefixes: &[&str], enabled: bool) -> io::Result<Self> {
        let mut pending = HashMap::new();
        if enabled {
            let mut available = available_memory();
            let mut jobs = Vec::new();
            for prefix in prefixes {
                let path = index_path(prefix);
                let bytes = fs::metadata(&path)?.len();
                let job = match available {
                    Some(free) if bytes > free => Pending::Done(Preload::Skipped {
                        path,
                        bytes,
                        available: free,
                    }),
                    _ => {
                        if let Some(free) = available.as_mut() {
                            *free -= bytes;
                        }
                        let (sender, receiver) = mpsc::channel();
                        jobs.push((path, sender));
                        Pending::Running(receiver)
                    }
                };
                pending.insert((*prefix).to_owned(), job);
            }
            if !jobs.is_empty() {
                thread::spawn(move || {
                    for (path, sender) in jobs {
                        // A dropped receiver only means nobody waits for it.
                        let _ = sender.send(preload(path));
                    }
                });
            }
        }
        Ok(Self {
            pending: Mutex::new(pending),
        })
    }

    /// Wait for the warm-up of `prefix` and report it. Reports `Disabled`
    /// for a prefix that was not started or was already waited for.
    pub(crate) fn wait(&self, prefix: &str) -> io::Result<Preload> {
        let job = self.pending.lock().expect("warm-up lock").remove(prefix);
        match job {
            Some(Pending::Running(receiver)) => receiver
                .recv()
                .unwrap_or_else(|_| Err(io::Error::other("index warm-up thread stopped"))),
            Some(Pending::Done(preload)) => Ok(preload),
            None => Ok(Preload::Disabled),
        }
    }
}

fn preload(path: String) -> io::Result<Preload> {
    let started = Instant::now();
    let mut file = File::open(&path)?;
    let mut buffer = vec![0; 16 << 20];
    let mut bytes = 0;
    loop {
        match file.read(&mut buffer)? {
            0 => break,
            read => bytes += read as u64,
        }
    }
    Ok(Preload::Loaded {
        path,
        bytes,
        seconds: started.elapsed().as_secs_f64(),
    })
}

/// Memory this process can fill with page cache: the kernel's
/// `MemAvailable`, capped by the headroom of every enclosing cgroup (Slurm
/// jobs, containers). File cache already in a cgroup counts as free, since
/// the kernel reclaims it first. `None` where this is not known (no `/proc`,
/// as on macOS).
fn available_memory() -> Option<u64> {
    available_memory_under(Path::new("/"))
}

/// `available_memory` with `/proc` and `/sys` under `root`.
fn available_memory_under(root: &Path) -> Option<u64> {
    let meminfo = fs::read_to_string(root.join("proc/meminfo")).ok()?;
    let mem_available = meminfo
        .lines()
        .find_map(|line| line.strip_prefix("MemAvailable:"))
        .and_then(|value| {
            value
                .trim()
                .trim_end_matches("kB")
                .trim()
                .parse::<u64>()
                .ok()
        })
        .map(|kilobytes| kilobytes * 1024);
    cgroup_headroom(root).into_iter().chain(mem_available).min()
}

/// Headroom of this process's memory cgroup and its ancestors, which hold
/// the limits Slurm sets on a job.
fn cgroup_headroom(root: &Path) -> Vec<u64> {
    let read_u64 = |path: PathBuf| fs::read_to_string(path).ok()?.trim().parse::<u64>().ok();
    let stat = |dir: &Path, key: &str| {
        fs::read_to_string(dir.join("memory.stat"))
            .ok()?
            .lines()
            .find_map(|line| {
                line.strip_prefix(key)?
                    .strip_prefix(' ')?
                    .trim()
                    .parse::<u64>()
                    .ok()
            })
    };
    let Ok(cgroups) = fs::read_to_string(root.join("proc/self/cgroup")) else {
        return Vec::new();
    };
    let mut headroom = Vec::new();
    for line in cgroups.lines() {
        let mut fields = line.splitn(3, ':');
        let (Some(_), Some(controllers), Some(path)) =
            (fields.next(), fields.next(), fields.next())
        else {
            continue;
        };
        // cgroup v2 has one unified hierarchy; v1 has a "memory" controller.
        let (mount, limit, usage, cache) = if controllers.is_empty() {
            ("sys/fs/cgroup", "memory.max", "memory.current", "file")
        } else if controllers
            .split(',')
            .any(|controller| controller == "memory")
        {
            (
                "sys/fs/cgroup/memory",
                "memory.limit_in_bytes",
                "memory.usage_in_bytes",
                "total_cache",
            )
        } else {
            continue;
        };
        let mount = root.join(mount);
        let mut dir = Some(mount.join(path.trim_start_matches('/')));
        while let Some(current) = dir {
            // "max" (v2) does not parse; v1 reports no limit as a huge value.
            if let (Some(limit), Some(usage)) = (
                read_u64(current.join(limit)).filter(|limit| *limit < 1 << 60),
                read_u64(current.join(usage)),
            ) {
                let reclaimable = stat(&current, cache).unwrap_or(0);
                headroom.push(limit.saturating_sub(usage.saturating_sub(reclaimable)));
            }
            dir = current
                .parent()
                .filter(|parent| parent.starts_with(&mount))
                .map(Path::to_path_buf);
        }
    }
    headroom
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn warmup_reads_or_skips_indexes() {
        let dir = std::env::temp_dir().join(format!("aeon-warmup-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let prefix = dir.join("nt").to_string_lossy().into_owned();
        let second = dir.join("aa").to_string_lossy().into_owned();
        fs::write(index_path(&prefix), vec![7u8; 1 << 20]).unwrap();
        fs::write(index_path(&second), vec![7u8; 1 << 10]).unwrap();

        let warmup = Warmup::start(&[&prefix, &second], true).unwrap();
        // Small files fit anywhere these tests run; wait out of order.
        for (name, size) in [(&second, 1 << 10), (&prefix, 1 << 20)] {
            match warmup.wait(name).unwrap() {
                Preload::Loaded { bytes, .. } => assert_eq!(bytes, size),
                other => panic!("unexpected {other:?}"),
            }
        }
        // A finished or unknown warm-up reports nothing to wait for.
        assert_eq!(warmup.wait(&prefix).unwrap(), Preload::Disabled);
        let off = Warmup::start(&[&prefix], false).unwrap();
        assert_eq!(off.wait(&prefix).unwrap(), Preload::Disabled);
        fs::remove_dir_all(&dir).unwrap();
    }

    fn write(root: &Path, path: &str, contents: &str) {
        let path = root.join(path);
        fs::create_dir_all(path.parent().unwrap()).unwrap();
        fs::write(path, contents).unwrap();
    }

    #[test]
    fn available_memory_takes_tightest_cgroup_limit() {
        let gib = 1u64 << 30;
        let base = std::env::temp_dir().join(format!("aeon-memory-{}", std::process::id()));
        let meminfo = format!(
            "MemTotal: {} kB\nMemAvailable: {} kB\n",
            64 * gib / 1024,
            40 * gib / 1024
        );

        // cgroup v2 (Slurm job and step): the job's 16 GiB limit, with 10 GiB
        // used of which 4 GiB is file cache, leaves 10 GiB; the step has no limit.
        let v2 = base.join("v2");
        write(&v2, "proc/meminfo", &meminfo);
        write(&v2, "proc/self/cgroup", "0::/job_1/step_0\n");
        write(
            &v2,
            "sys/fs/cgroup/job_1/memory.max",
            &format!("{}\n", 16 * gib),
        );
        write(
            &v2,
            "sys/fs/cgroup/job_1/memory.current",
            &format!("{}\n", 10 * gib),
        );
        write(
            &v2,
            "sys/fs/cgroup/job_1/memory.stat",
            &format!("anon 1\nfile {}\n", 4 * gib),
        );
        write(&v2, "sys/fs/cgroup/job_1/step_0/memory.max", "max\n");
        write(&v2, "sys/fs/cgroup/job_1/step_0/memory.current", "5\n");
        assert_eq!(available_memory_under(&v2), Some(10 * gib));

        // cgroup v1 without a limit falls back to MemAvailable.
        let v1 = base.join("v1");
        write(&v1, "proc/meminfo", &meminfo);
        write(
            &v1,
            "proc/self/cgroup",
            "5:cpuacct,cpu:/\n4:memory:/slurm/uid_1\n",
        );
        write(
            &v1,
            "sys/fs/cgroup/memory/slurm/uid_1/memory.limit_in_bytes",
            "9223372036854771712\n",
        );
        write(
            &v1,
            "sys/fs/cgroup/memory/slurm/uid_1/memory.usage_in_bytes",
            "100\n",
        );
        assert_eq!(available_memory_under(&v1), Some(40 * gib));

        // No /proc (macOS): unknown, so every index is read.
        assert_eq!(available_memory_under(&base.join("none")), None);
        fs::remove_dir_all(&base).unwrap();
    }
}
