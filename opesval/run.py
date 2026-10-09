"""Prepare and execute the simulations, one directory per replica.

Layout: runs/<set>/<tier>/<system>/<method>/<replica>/
where <set> defaults to the PLUMED install label.
"""

import json
import os
import re
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .analysis import analyze_run
from .systems import ROOT

OPT = Path(os.environ.get("OPT_DIR", ROOT / "opt"))
RUNS = Path(os.environ.get("RUNS_DIR", ROOT / "runs"))
DEFAULT_GROMACS = "2024.3"


def read_build_info(prefix):
    info = {}
    with open(prefix / "BUILD_INFO") as f:
        for line in f:
            key, _, value = line.strip().partition("=")
            info[key] = value
    return info


def plumed_install(label):
    prefix = OPT / "plumed" / label
    if not (prefix / "bin" / "plumed").exists():
        sys.exit(f"PLUMED '{label}' not installed in {prefix}, use build/build_plumed.sh")
    info = read_build_info(prefix)
    info["prefix"] = str(prefix)
    return info


def gromacs_install(version):
    prefix = OPT / "gromacs" / version
    if not (prefix / "bin" / "gmx").exists():
        sys.exit(f"GROMACS {version} not installed in {prefix}, use build/build_gromacs.sh")
    info = read_build_info(prefix)
    info["prefix"] = str(prefix)
    return info


OPES_ACTION = re.compile(r"^\s*\w+:\s*OPES_METAD(_EXPLORE)?\b")


def add_opes_keywords(text, extra):
    """Append keywords, e.g. "CUTOFF_SCHEME=PR", to the OPES action of a PLUMED input."""
    if not extra:
        return text
    lines = text.splitlines(keepends=True)
    for n, line in enumerate(lines):
        if OPES_ACTION.match(line):
            if line.rstrip().endswith("..."):  # multi-line action
                lines.insert(n + 1, f"  {extra}\n")
            else:
                lines[n] = line.rstrip("\n") + f" {extra}\n"
            return "".join(lines)
    raise ValueError("no OPES_METAD or OPES_METAD_EXPLORE action found")


def default_run_set(label, opes_extra):
    if not opes_extra:
        return label
    return label + "+" + "+".join(kw.lower().replace("=", "-") for kw in opes_extra.split())


@dataclass
class Job:
    system: object
    method: object
    tier: str
    replica: int
    run_set: str
    plumed: dict
    gromacs: dict = None
    opes_extra: str = ""  # keywords added to the OPES action

    @property
    def nsteps(self):
        return self.system.tiers[self.tier].nsteps

    @property
    def dir(self):
        return RUNS / self.run_set / self.tier / self.system.name / self.method.name / f"{self.replica:03d}"

    @property
    def name(self):
        return f"{self.run_set}/{self.tier}/{self.system.name}/{self.method.name}/{self.replica:03d}"

    def signature(self):
        """What must match for an existing run to be reused."""
        sig = {"plumed_sha": self.plumed["sha"], "nsteps": self.nsteps}
        sig["plumed_input"] = self.plumed_input()
        if self.gromacs:
            sig["gromacs"] = self.gromacs["version"]
        return sig

    def plumed_input(self):
        return add_opes_keywords((self.system.path / self.method.plumed).read_text(), self.opes_extra)

    def is_done(self):
        try:
            with open(self.dir / "run.json") as f:
                info = json.load(f)
        except FileNotFoundError:
            return False
        return info.get("status") == "done" and info.get("signature") == self.signature()


def make_jobs(pairs, tier, replicas, labels, run_set=None, gromacs_version=DEFAULT_GROMACS, opes_extra=""):
    jobs = []
    gromacs = None
    plumeds = {label: plumed_install(label) for label in labels}
    for system, method in pairs:
        if tier not in system.tiers:
            sys.exit(f"tier '{tier}' not defined for {system.name}")
        if system.engine == "gromacs" and gromacs is None:
            gromacs = gromacs_install(gromacs_version)
        nrep = system.tiers[tier].replicas
        for replica in replicas if replicas is not None else range(nrep):
            if replica >= available_replicas(system):
                sys.exit(f"replica {replica} not available for {system.name}")
            for label in labels:  # interleave the PLUMED versions
                jobs.append(
                    Job(system, method, tier, replica, run_set or default_run_set(label, opes_extra), plumeds[label],
                        gromacs if system.engine == "gromacs" else None, opes_extra)
                )
    # longest first, for a better load balance
    jobs.sort(key=lambda j: -j.nsteps)
    return jobs


def available_replicas(system):
    if system.engine == "ves":
        return len(np.loadtxt(system.path / "replicas.dat", ndmin=2))
    return len(list((system.path / "tpr").glob("input.*.tpr")))


def prepare(job):
    if job.dir.exists():
        shutil.rmtree(job.dir)
    job.dir.mkdir(parents=True)
    (job.dir / "plumed.dat").write_text(job.plumed_input())
    env = dict(os.environ, OMP_NUM_THREADS="1", PLUMED_NUM_THREADS="1")
    if job.system.engine == "ves":
        lines = (job.system.path / "replicas.dat").read_text().splitlines()
        rows = [line.split() for line in lines if not line.startswith("#")]
        _, seed, x0, y0 = rows[job.replica]  # keep the values verbatim
        template = (job.system.path / "md_input.dat").read_text()
        for key, value in {"@NSTEP@": job.nsteps, "@SEED@": seed, "@X0@": x0, "@Y0@": y0}.items():
            template = template.replace(key, str(value))
        (job.dir / "md_input.dat").write_text(template)
        shutil.copy(job.system.path / "potential.dat", job.dir / "potential.dat")
        cmd = [f"{job.plumed['prefix']}/bin/plumed", "ves_md_linearexpansion", "md_input.dat"]
        env["LD_LIBRARY_PATH"] = f"{job.plumed['prefix']}/lib:" + env.get("LD_LIBRARY_PATH", "")
    elif job.system.engine == "gromacs":
        shutil.copy(job.system.path / "tpr" / f"input.{job.replica}.tpr", job.dir / "input.tpr")
        cmd = [f"{job.gromacs['prefix']}/bin/gmx", "mdrun", "-s", "input.tpr", "-plumed", "plumed.dat",
               "-deffnm", "md", "-nsteps", str(job.nsteps), "-ntmpi", "1", "-ntomp", "1", "-pin", "off", "-nb", "cpu"]
        # GROMACS is patched in runtime mode, the kernel is chosen here
        env["PLUMED_KERNEL"] = f"{job.plumed['prefix']}/lib/libplumedKernel.so"
    else:
        raise ValueError(job.system.engine)
    return cmd, env


def check_plumed_version(job):
    """Make sure the run used the intended PLUMED kernel."""
    log = job.dir / ("plumed.log" if job.system.engine == "ves" else "md.log")
    match = re.search(r"Version: \S+ \(git: ([0-9a-f]+)\)", log.read_text(errors="replace"))
    if not match:
        raise RuntimeError(f"PLUMED version not found in {log}")
    if not job.plumed["sha"].startswith(match.group(1)):
        raise RuntimeError(f"wrong PLUMED kernel: git {match.group(1)}, expected {job.plumed['sha']}")
    return match.group(1)


def execute(job):
    info = {"job": job.name, "signature": job.signature(), "status": "running",
            "plumed": {k: job.plumed[k] for k in ("label", "ref", "sha")},
            "gromacs": job.gromacs and job.gromacs["version"], "start": time.strftime("%Y-%m-%dT%H:%M:%S")}
    t0 = time.time()
    try:
        cmd, env = prepare(job)
        info["command"] = cmd
        with open(job.dir / "out.log", "w") as out:
            proc = subprocess.run(cmd, cwd=job.dir, env=env, stdout=out, stderr=subprocess.STDOUT)
        info["wall_time"] = time.time() - t0
        if proc.returncode != 0:
            raise RuntimeError(f"exit code {proc.returncode}, see {job.dir / 'out.log'}")
        info["plumed_git"] = check_plumed_version(job)
        analyze_run(job.system, job.dir)
        info["status"] = "done"
    except Exception as e:  # keep going with the other jobs
        info["status"] = "failed"
        info["error"] = str(e)
        info.setdefault("wall_time", time.time() - t0)
    job.dir.mkdir(parents=True, exist_ok=True)
    with open(job.dir / "run.json", "w") as f:
        json.dump(info, f, indent=1)
    return job, info


def run_jobs(jobs, njobs, force=False, dry_run=False):
    todo = [j for j in jobs if force or not j.is_done()]
    print(f"{len(jobs)} jobs selected, {len(jobs) - len(todo)} already done, {len(todo)} to run on {njobs} cores")
    if dry_run:
        for j in todo:
            print(f"  {j.name}  nsteps={j.nsteps}")
        return True
    failed = 0
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=njobs) as pool:
        futures = [pool.submit(execute, j) for j in todo]
        for n, future in enumerate(as_completed(futures), 1):
            job, info = future.result()
            failed += info["status"] != "done"
            msg = info.get("error", "")
            elapsed = (time.time() - t0) / 60
            print(f"[{n}/{len(todo)} {elapsed:.0f}min] {job.name} {info['status']} ({info['wall_time']:.0f}s) {msg}",
                  flush=True)
    if failed:
        print(f"{failed} jobs failed")
    return failed == 0


def reanalyze(jobs):
    for job in jobs:
        if job.is_done():
            analyze_run(job.system, job.dir)
            print(f"analyzed {job.name}")
