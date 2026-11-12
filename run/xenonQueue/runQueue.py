#!/usr/bin/env python3
"""Run a queue of xenonTransferInterFoam cases one after the other and record each run.

Usage (OpenFOAM and repository environment sourced: <openfoam>/etc/bashrc, <repo>/etc/sourceme.sh):
    python3 runQueue.py <queue.py> --np N [--only PATTERN ...] [--force] [--prepare-only]
    python3 runQueue.py <queue.py> --list
    sbatch -n N [-t HH:MM:SS] runQueue.slurm <queue.py> [--only PATTERN ...]   (from this folder)

The queue file defines RUNS, a list of dicts:
    name   run folder under RUNS_DIR (default: runs/ next to the queue file)
    case   template case under tom-v2606/run (xenonDroplet, xenonDropStatic)
    args   Allrun options, e.g. ["-init", "smooth", "-cubic"]
    set    {file: {keyword: value}}: every line "keyword value;" of the file gets the value
    sub    [(file, regex, replacement)]: any other edit (re.M)
    note   free text, copied to the record
Every edit must match at least once, or the run fails before it starts.

Each run: copy the template, apply the edits (plus numberOfSubdomains = N), Allrun -prepare (mesh, fields,
decomposePar), the solver (srun under SLURM, mpirun otherwise; log appended), reconstructPar. The record is
<run>/runInfo.json and one row of RUNS_DIR/summary.csv (rebuilt from all records after each run).
Finished runs are skipped (--force reruns them). A run cut by the job time limit is resumed from its latest
time (startFrom latestTime) when the queue is submitted again; the boundary outflow (XEOUT, static case) then
counts the last segment only. A failed run is recorded and the queue goes on.
"""
import argparse
import csv
import datetime
import fnmatch
import json
import math
import os
import re
import runpy
import shutil
import socket
import subprocess
import sys
import time
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "src" / "xenonTransfer" / "test" / "newmanSphere"))
import checkNewmanCase as cn          # noqa: E402  Newman fraction remaining

SKIP = ("processor", "log.", "postProcessing", "dynamicCode", "analysis", "runInfo.json")


def now():
    return datetime.datetime.now().isoformat(timespec="seconds")


def is_time(name):
    try:
        float(name)
        return True
    except ValueError:
        return False


def copy_template(src, dst):
    # case set-up only: no results, logs or decomposition
    def ignore(_, names):
        return [n for n in names if n.startswith(SKIP) or is_time(n)]
    shutil.copytree(src, dst, ignore=ignore)


def apply_edits(d, R, np_):
    done = []
    sets = {f: dict(kv) for f, kv in R.get("set", {}).items()}
    sets.setdefault("system/decomposeParDict", {})["numberOfSubdomains"] = np_
    for f, kv in sets.items():
        p = d / f
        s = p.read_text()
        for k, v in kv.items():
            s, n = re.subn(r"^(\s*)%s(\s+)[^;\n]*;" % re.escape(k), r"\g<1>%s\g<2>%s;" % (k, v), s, flags=re.M)
            if n == 0:
                raise ValueError("%s: no entry '%s'" % (f, k))
            done.append("%s: %s %s (%d)" % (f, k, v, n))
        p.write_text(s)
    for f, pat, rep in R.get("sub", []):
        p = d / f
        s, n = re.subn(pat, rep, p.read_text(), flags=re.M)
        if n == 0:
            raise ValueError("%s: no match for %s" % (f, pat))
        p.write_text(s)
        done.append("%s: %s -> %s (%d)" % (f, pat, rep, n))
    return done


def sh(cmd, cwd):
    return subprocess.run(["bash", "-c", cmd], cwd=str(cwd)).returncode


def entry(path, key):
    m = re.search(r"^\s*%s\s+([^;\s]+)\s*;" % re.escape(key), re.sub(r"//.*", "", path.read_text()), re.M)
    return m.group(1) if m else None


# ---------------------------------------------------------------------------- results

def series(d, name):
    rows = {}
    for p in sorted((d / "postProcessing" / name).glob("*/volFieldValue*.dat"),
                    key=lambda q: float(q.parent.name)):
        for line in open(p):
            if line.strip() and not line.startswith("#"):
                a = line.split()
                rows[float(a[0])] = float(a[1])
    return sorted(rows.items())


def umax(d):
    rows = {}
    for p in sorted((d / "postProcessing" / "UminMax").glob("*/fieldMinMax*.dat"),
                    key=lambda q: float(q.parent.name)):
        for line in open(p):
            if line.strip() and not line.startswith("#"):
                a = re.sub(r"\([^)]*\)", " ", line).split()
                rows[float(a[0])] = float(a[4])
    return sorted(rows.items())


def parse_log(log):
    # steps, end time, wall time (summed over restarts), drop 0, XEOUT, capacity limiter, fatal error
    L = dict(steps=0, t=None, dtMax=0.0, wall=0.0, cap=0, capRejectedMax=0.0, fatal=False, end=False)
    drop0, xe, segWall, nDropsMax, seen = [], None, 0.0, 0, set()
    t, dt, restart = None, 0.0, False
    if not log.is_file():
        return L, drop0
    for line in open(log, errors="replace"):
        if line.startswith("Time = "):
            t = float(line.split()[2])
            if restart:                           # resumed: drop the steps after the restart time
                seen = {s for s in seen if s < t - 0.5 * dt}
                restart = False
            seen.add(t)
        elif line.startswith("deltaT = "):
            dt = float(line.split()[2])
            L["dtMax"] = max(L["dtMax"], dt)
        elif line.startswith("Starting time loop"):
            restart = t is not None
        elif line.startswith("ExecutionTime = "):
            w = float(line.split()[2])
            if w < segWall:                       # restarted: new segment
                L["wall"] += segWall
            segWall = w
        elif line.startswith("  drop 0:") and t is not None:
            g = dict(re.findall(r" (\w+)=([-+0-9.eE]+)(?=\s)", line))
            rec = dict(t=t, **{k: float(g[k]) for k in ("age", "Dvol", "Ufall", "AR") if k in g})
            m = re.search(r" Ccm=\(([^)]*)\)", line)
            if m:
                rec["z"] = float(m.group(1).split()[2])
            if drop0 and drop0[-1]["t"] == t:
                drop0[-1] = rec
            else:
                drop0.append(rec)
        elif line.startswith("gasDiffusionSimple:"):
            m = re.search(r" nDrops=(\d+)", line)
            if m:
                nDropsMax = max(nDropsMax, int(m.group(1)))
        elif line.startswith("XEOUT "):
            xe = float(line.split()[3])
        elif line.startswith("MT CAP SUMMARY"):
            L["cap"] += 1
            m = re.search(r"rejectedAbsRatio=([-+0-9.eE]+)", line)
            if m:
                L["capRejectedMax"] = max(L["capRejectedMax"], float(m.group(1)))
        elif "FOAM FATAL" in line or "Floating point exception" in line and "trapping" not in line:
            L["fatal"] = True
        elif line.strip() == "End":
            L["end"] = True
    L["t"], L["wall"], L["xeOut"], L["nDropsMax"] = t, L["wall"] + segWall, xe, nDropsMax
    L["steps"] = len(seen)
    return L, drop0


def collect(d):
    app = entry(d / "system" / "controlDict", "application")
    L, drop0 = parse_log(d / ("log." + app))
    r = dict(endTimeReached=L["t"], steps=L["steps"], dtMax=L["dtMax"], solverWall_h=L["wall"] / 3600,
             capEvents=L["cap"], capRejectedMax=L["capRejectedMax"], nDropsMax=L["nDropsMax"])
    pp = d / "constant" / "phaseProperties"
    Dliq = float(entry(pp, "Dliq") or "nan")
    r["Dliq"] = Dliq
    cm = d / "log.checkMesh"
    if cm.is_file():
        s = cm.read_text()
        m = re.search(r"^\s*cells:\s+(\d+)", s, re.M)
        r["cells"] = int(m.group(1)) if m else None
        m = re.search(r"Min volume = ([-+0-9.eE]*[0-9])", s)
        r["dxMin_um"] = 1e6 * float(m.group(1)) ** (1 / 3) if m else None
    liq, tot = series(d, "ZXeMass_in_liquid"), dict(series(d, "ZXeMass_total"))
    if liq and tot:
        M0 = tot[min(tot)]
        tEnd, MlEnd = liq[-1]
        r.update(M0=M0, F_end=1 - MlEnd / M0, Mtot_end_M0=tot[max(tot)] / M0)
        if L["xeOut"] is not None:
            r["xeOut_M0"] = L["xeOut"] / M0
        if drop0 and "age" in drop0[0] and math.isfinite(Dliq):
            d0 = drop0[0]["Dvol"]
            ts, ages = [q["t"] for q in drop0], [q["age"] for q in drop0]

            def age(t):
                if t <= ts[0]:
                    return ages[0] * t / ts[0] if ts[0] > 0 else 0.0
                i = next((k for k in range(1, len(ts)) if ts[k] >= t), len(ts) - 1)
                f = (t - ts[i - 1]) / max(ts[i] - ts[i - 1], 1e-300)
                return ages[i - 1] + f * (ages[i] - ages[i - 1])
            err = [Ml / M0 / cn.remaining(4 * Dliq * age(t) / d0 ** 2) - 1 for t, Ml in liq]
            r.update(d0_um=1e6 * d0, theta_end=age(tEnd), Fo_end=4 * Dliq * age(tEnd) / d0 ** 2,
                     F_newman_end=1 - cn.remaining(4 * Dliq * age(tEnd) / d0 ** 2),
                     newmanErrMax=max(err, key=abs), newmanErrEnd=err[-1])
    if drop0:
        late = [q for q in drop0 if q["t"] >= 0.8 * drop0[-1]["t"]]
        if "Ufall" in late[0]:
            r["Ufall_late"] = sum(q["Ufall"] for q in late) / len(late)
        if "z" in drop0[-1]:
            r["zDrop_end_mm"] = 1e3 * drop0[-1]["z"]
        if "AR" in drop0[-1]:
            r["AR_end"] = drop0[-1]["AR"]
    U = umax(d)
    if U:
        late = [u for t, u in U if t >= U[0][0] + 0.9 * (U[-1][0] - U[0][0])]
        r.update(Umax_peak=max(u for _, u in U), Umax_late=sum(late) / len(late))
    return r, L


# ---------------------------------------------------------------------------- queue

def record(d, info):
    (d / "runInfo.json").write_text(json.dumps(info, indent=1, default=str))


def write_summary(runs_dir, names):
    rows = []
    for n in names:
        p = runs_dir / n / "runInfo.json"
        if p.is_file():
            info = json.loads(p.read_text())
            row = {k: v for k, v in info.items() if not isinstance(v, (dict, list))}
            row["args"] = " ".join(info.get("args", []))
            row.update(info.get("results", {}))
            rows.append(row)
    if not rows:
        return
    cols = []
    for row in rows:
        cols += [k for k in row if k not in cols]
    with open(runs_dir / "summary.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)


def git_version():
    try:
        h = subprocess.run(["git", "-C", str(REPO), "rev-parse", "--short", "HEAD"], capture_output=True,
                           text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", str(REPO), "status", "--porcelain", "--", "src/xenonTransfer",
                                "applications/solvers/xenonTransferInterFoam"], capture_output=True,
                               text=True).stdout.strip()
        return h + ("+dirty" if dirty else "")
    except OSError:
        return "?"


def run_one(R, runs_dir, np_, force, prepare_only):
    d = runs_dir / R["name"]
    info_path = d / "runInfo.json"
    old = json.loads(info_path.read_text()) if info_path.is_file() else {}
    if old.get("status") == "done" and not force:
        print("[%s] done, skipped" % R["name"])
        return
    resume = (not force and old.get("status") == "running" and (d / "processor0").is_dir()
              and old.get("np") == np_)
    info = dict(name=R["name"], case=R["case"], args=R.get("args", []), note=R.get("note", ""), np=np_,
                host=socket.gethostname(), job=os.environ.get("SLURM_JOB_ID", ""), version=git_version(),
                started=old.get("started") if resume else now(), segments=old.get("segments", 0) + 1 if resume else 1,
                edits=old.get("edits", []), status="preparing")
    t0 = time.time()
    try:
        if not resume:
            if d.exists():
                shutil.rmtree(d)
            copy_template(REPO / "run" / R["case"], d)
            info["edits"] = apply_edits(d, R, np_)
            record(d, info)
            print("[%s] prepare: Allrun -prepare %s" % (R["name"], " ".join(R.get("args", []))), flush=True)
            if sh("./Allrun -prepare " + " ".join(R.get("args", [])), d):
                raise RuntimeError("Allrun -prepare failed (see the log.* files)")
        if prepare_only:
            info["status"] = "prepared"
            return
        app = entry(d / "system" / "controlDict", "application")
        launch = ("srun -n %d" % np_) if os.environ.get("SLURM_JOB_ID") else ("mpirun -np %d" % np_)
        launch = os.environ.get("XENON_MPIRUN", launch)
        info["status"] = "running"
        record(d, info)
        print("[%s] %s: %s %s -parallel" % (R["name"], "resume" if resume else "run", launch, app), flush=True)
        rc = sh("%s %s -parallel >> log.%s 2>&1" % (launch, app, app), d)
        sh("reconstructPar -newTimes > log.reconstructPar 2>&1", d)
        res, L = collect(d)
        info["results"] = res
        endTime = float(entry(d / "system" / "controlDict", "endTime"))
        reached = L["t"] is not None and L["t"] >= endTime * (1 - 1e-6) - 1e-12
        info["status"] = "done" if rc == 0 and L["end"] and reached else "failed"
        info["solverRc"] = rc
        if info["status"] == "failed":
            log = d / ("log." + app)
            info["logTail"] = open(log, errors="replace").read().splitlines()[-25:] if log.is_file() else []
    except Exception as e:                                   # noqa: BLE001  record and go on with the queue
        info["status"], info["error"] = "failed", str(e)
    finally:
        info["finished"] = now()
        info["segmentWall_h"] = (time.time() - t0) / 3600
        if d.is_dir():
            record(d, info)
        print("[%s] %s (%.2f h)%s" % (R["name"], info["status"], info["segmentWall_h"],
                                      ": " + info["error"] if "error" in info else ""), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("queue")
    ap.add_argument("--np", type=int, default=int(os.environ.get("SLURM_NTASKS", "0") or 0))
    ap.add_argument("--only", nargs="*", default=None, help="run names or glob patterns")
    ap.add_argument("--force", action="store_true", help="rerun finished runs from scratch")
    ap.add_argument("--prepare-only", action="store_true", help="set up the runs (mesh, fields), no solver")
    ap.add_argument("--list", action="store_true", help="show the queue and the state of each run")
    a = ap.parse_args()

    qfile = Path(a.queue).resolve()
    Q = runpy.run_path(str(qfile))
    runs_dir = Path(Q.get("RUNS_DIR", qfile.parent / "runs"))
    runs = [R for R in Q["RUNS"] if a.only is None or any(fnmatch.fnmatch(R["name"], p) for p in a.only)]
    if len({R["name"] for R in Q["RUNS"]}) != len(Q["RUNS"]):
        sys.exit("run names in %s are not unique" % qfile)

    if a.list:
        for R in runs:
            p = runs_dir / R["name"] / "runInfo.json"
            info = json.loads(p.read_text()) if p.is_file() else {}
            res = info.get("results", {})
            print("%-34s %-10s %6s h  t = %-10s %s" % (R["name"], info.get("status", "-"),
                  "%.2f" % res["solverWall_h"] if "solverWall_h" in res else "-",
                  res.get("endTimeReached", "-"), R.get("note", "")))
        return
    if a.np < 2:
        sys.exit("give --np N (N >= 2), or submit with sbatch -n N")
    if os.environ.get("FOAM_USER_APPBIN"):
        os.environ["PATH"] = os.environ["FOAM_USER_APPBIN"] + os.pathsep + os.environ["PATH"]
    if "WM_PROJECT_DIR" not in os.environ or shutil.which("xenonTransferInterFoam") is None:
        sys.exit("source the OpenFOAM bashrc and <repo>/etc/sourceme.sh first (xenonTransferInterFoam not found)")

    runs_dir.mkdir(parents=True, exist_ok=True)
    print("%s: %d runs on %d processes, results in %s" % (qfile.name, len(runs), a.np, runs_dir), flush=True)
    for R in runs:
        run_one(R, runs_dir, a.np, a.force, a.prepare_only)
        write_summary(runs_dir, [q["name"] for q in Q["RUNS"]])
    print("summary: %s" % (runs_dir / "summary.csv"))


if __name__ == "__main__":
    main()
