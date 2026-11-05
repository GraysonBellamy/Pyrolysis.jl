"""Shared measurement protocol for the timing scripts in this directory.

Every timed process runs pinned to one performance core (CPU 1 by default; its
hyperthread sibling must stay idle), one at a time, after a quiet-machine gate.
Each run records what could have disturbed it: the busy fraction of the timed
core, of its sibling and of every other core over the run (from /proc/stat),
the involuntary context switches of the process, the one-minute load and the
package temperature before it starts. A run is judged afterwards on those
records, never discarded silently.

Python 3.11+, standard library only.
"""

from pathlib import Path
import hashlib
import json
import os
import platform
import subprocess
import threading
import time

CPU = int(os.environ.get("TIMING_CPU", "1"))


def keep_driver_off(cpu=CPU):
    """Move the calling (driver) process off the timed core and its siblings,
    onto the efficiency cores when there are any."""
    avoid = {cpu, *sibling(cpu)}
    pcores = performance_cores()
    allowed = set(range(os.cpu_count())) - avoid
    if pcores:
        e_cores = allowed - pcores
        allowed = e_cores or allowed
    os.sched_setaffinity(0, allowed)
    return sorted(allowed)


def _read(path, default=None):
    try:
        return Path(path).read_text().strip()
    except OSError:
        return default


def sibling(cpu=CPU):
    """Hyperthread siblings of `cpu`, excluding it."""
    text = _read(f"/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list", str(cpu))
    out = set()
    for part in text.split(","):
        a, _, b = part.partition("-")
        out.update(range(int(a), int(b or a) + 1))
    out.discard(cpu)
    return sorted(out)


def performance_cores():
    text = _read("/sys/devices/cpu_core/cpus")
    if not text:
        return None
    out = set()
    for part in text.split(","):
        a, _, b = part.partition("-")
        out.update(range(int(a), int(b or a) + 1))
    return out


def cpu_times():
    """{cpu: (busy_jiffies, total_jiffies)} from /proc/stat."""
    out = {}
    for line in Path("/proc/stat").read_text().splitlines():
        if not line.startswith("cpu") or line.startswith("cpu "):
            continue
        name, *fields = line.split()
        v = [int(x) for x in fields]
        idle = v[3] + v[4]  # idle + iowait
        out[int(name[3:])] = (sum(v[:8]) - idle, sum(v[:8]))
    return out


def busy_fractions(before, after, cpu=CPU):
    """Busy fraction of the timed core, its sibling(s) and the mean of the rest."""

    def frac(c):
        db = after[c][0] - before[c][0]
        dt = after[c][1] - before[c][1]
        return db / dt if dt > 0 else 0.0

    sib = sibling(cpu)
    rest = [c for c in after if c != cpu and c not in sib]
    return dict(
        cpu_busy=frac(cpu),
        sibling_busy=max((frac(c) for c in sib), default=0.0),
        others_busy=sum(frac(c) for c in rest) / max(len(rest), 1),
    )


def package_temperature():
    for hw in Path("/sys/class/hwmon").glob("hwmon*"):
        if _read(hw / "name") == "coretemp":
            t = _read(hw / "temp1_input")
            return int(t) / 1000 if t else None
    return None


def wait_quiet(max_load=2.0, max_sibling=0.05, max_others=0.10, timeout=1800, log=print):
    """Block until the one-minute load is below `max_load` and, over 0.3 s,
    the timed core's sibling and the other cores are idle enough."""
    start = time.monotonic()
    while True:
        before = cpu_times()
        time.sleep(0.3)
        b = busy_fractions(before, cpu_times())
        load = os.getloadavg()[0]
        if load < max_load and b["sibling_busy"] <= max_sibling and b["others_busy"] <= max_others:
            return
        if time.monotonic() - start > timeout:
            raise RuntimeError(f"machine not quiet after {timeout} s: load {load}, {b}")
        log(f"waiting for a quiet machine: load {load:.2f}, sibling {b['sibling_busy']:.2f}, "
            f"others {b['others_busy']:.2f}")
        time.sleep(4.0)


def run_pinned(command, cwd, log_path, env=None, stdin_text=None, timeout=None, cpu=CPU):
    """Run `command` pinned to `cpu`; return its wall time, resource usage and
    the disturbance record. `taskset` execs the command, so the pid is the
    program's own and wait4's rusage is the program's."""
    load_before = os.getloadavg()
    temp_before = package_temperature()
    with open(log_path, "w") as log:
        t0 = time.perf_counter_ns()
        c0 = cpu_times()
        proc = subprocess.Popen(
            ["taskset", "-c", str(cpu), *command],
            cwd=cwd,
            env=env,
            stdin=subprocess.PIPE if stdin_text is not None else subprocess.DEVNULL,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
        if stdin_text is not None:
            proc.stdin.write(stdin_text)
            proc.stdin.close()
        timer = None
        if timeout is not None:
            timer = threading.Timer(timeout, proc.kill)
            timer.start()
        _, status, ru = os.wait4(proc.pid, 0)  # blocking: the driver stays asleep
        t1 = time.perf_counter_ns()
        c1 = cpu_times()
        if timer is not None:
            timer.cancel()
    proc.returncode = os.waitstatus_to_exitcode(status)
    return dict(
        returncode=proc.returncode,
        timed_out=timeout is not None and proc.returncode == -9,
        wall_s=(t1 - t0) / 1e9,
        user_s=ru.ru_utime,
        sys_s=ru.ru_stime,
        maxrss_kib=ru.ru_maxrss,
        voluntary_switches=ru.ru_nvcsw,
        involuntary_switches=ru.ru_nivcsw,
        load_before=load_before[0],
        temperature_before_c=temp_before,
        **busy_fractions(c0, c1, cpu),
    )


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def environment(cpu=CPU):
    """What the machine was doing and how it was configured."""
    cpu_dir = Path(f"/sys/devices/system/cpu/cpu{cpu}/cpufreq")
    model = next(
        (x.split(":", 1)[1].strip() for x in Path("/proc/cpuinfo").read_text().splitlines()
         if x.startswith("model name")),
        platform.processor(),
    )
    ac = [_read(p) for p in Path("/sys/class/power_supply").glob("A*/online")]
    pcores = performance_cores()
    return dict(
        time=time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        host=platform.node(),
        kernel=platform.release(),
        cpu_model=model,
        logical_cpus=os.cpu_count(),
        timed_cpu=cpu,
        timed_cpu_sibling=sibling(cpu),
        timed_cpu_is_performance_core=None if pcores is None else cpu in pcores,
        max_frequency_khz=_read(cpu_dir / "cpuinfo_max_freq"),
        governor=_read(cpu_dir / "scaling_governor"),
        energy_performance_preference=_read(cpu_dir / "energy_performance_preference"),
        no_turbo=_read("/sys/devices/system/cpu/intel_pstate/no_turbo"),
        platform_profile=_read("/sys/firmware/acpi/platform_profile"),
        ac_online=ac,
        loadavg=os.getloadavg(),
        temperature_c=package_temperature(),
        uptime_s=float(_read("/proc/uptime", "0 0").split()[0]),
    )


def check_environment(env):
    """Refuse a configuration that would bias the times."""
    problems = []
    if env["timed_cpu_is_performance_core"] is False:
        problems.append(f"CPU {env['timed_cpu']} is not a performance core")
    if env["ac_online"] and "1" not in env["ac_online"]:
        problems.append("on battery")
    if env["no_turbo"] == "1":
        problems.append("turbo disabled")
    if env["platform_profile"] not in (None, "performance"):
        problems.append(f"platform profile {env['platform_profile']}")
    if problems:
        raise RuntimeError("timing environment: " + "; ".join(problems))


def git_stamp(repo):
    def g(*a):
        return subprocess.check_output(["git", "-C", str(repo), *a], text=True).strip()

    return dict(commit=g("rev-parse", "HEAD"), tree=g("rev-parse", "HEAD^{tree}"),
                dirty=bool(g("status", "--porcelain")))


def append_jsonl(path, row):
    with open(path, "a") as f:
        f.write(json.dumps(row) + "\n")


def read_jsonl(path):
    p = Path(path)
    if not p.exists():
        return []
    return [json.loads(x) for x in p.read_text().splitlines() if x.strip()]
