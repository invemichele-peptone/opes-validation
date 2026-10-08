"""Command line interface: python -m opesval {list,run,analyze,compare} ..."""

import argparse
import os
import sys

from . import systems
from .run import DEFAULT_GROMACS, RUNS, make_jobs, reanalyze, run_jobs


def parse_replicas(text):
    """'0-4,7' -> [0, 1, 2, 3, 4, 7]"""
    if text is None:
        return None
    out = []
    for part in text.split(","):
        lo, _, hi = part.partition("-")
        out += range(int(lo), int(hi or lo) + 1)
    return out


def split(text):
    return text.split(",") if text else None


def add_selection(p):
    p.add_argument("-s", "--systems", help="comma separated systems (default: all)")
    p.add_argument("-m", "--methods", help="comma separated methods (default: all non optional)")
    p.add_argument("-t", "--tier", default="quick", help="'quick' or 'full' (default: quick)")
    p.add_argument("--optional", action="store_true", help="include optional methods")


def cmd_list(args):
    for s in systems.SYSTEMS.values():
        tiers = ", ".join(f"{k}: {t.replicas}x{t.nsteps:.0e} steps" for k, t in s.tiers.items())
        print(f"{s.name} [{s.engine}] {s.description}\n  tiers  {tiers}")
        for m in s.methods:
            print(f"  method {m.name:10s} {m.plumed:22s} ({m.paper} paper){' optional' if m.optional else ''}")


def cmd_run(args):
    if args.set and len(args.plumed) > 1:
        sys.exit("--set can be used only with a single --plumed")
    pairs = list(systems.select(split(args.systems), split(args.methods), args.optional))
    jobs = make_jobs(pairs, args.tier, parse_replicas(args.replicas), args.plumed, args.set, args.gromacs)
    ok = run_jobs(jobs, args.jobs, force=args.force, dry_run=args.dry_run)
    sys.exit(0 if ok else 1)


def cmd_analyze(args):
    pairs = list(systems.select(split(args.systems), split(args.methods), args.optional))
    labels = args.plumed
    jobs = make_jobs(pairs, args.tier, parse_replicas(args.replicas), labels, args.set, args.gromacs)
    reanalyze(jobs)


def cmd_compare(args):
    from .compare import compare

    pairs = list(systems.select(split(args.systems), split(args.methods), True))
    outdir = args.out or os.path.join("reports", f"{args.A}_vs_{args.B}", args.tier)
    flagged = compare(args.A, args.B, args.tier, pairs, outdir)
    sys.exit(1 if flagged and args.strict else 0)


def main():
    parser = argparse.ArgumentParser(prog="opesval", description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("list", help="list systems, methods and tiers").set_defaults(func=cmd_list)

    p = sub.add_parser("run", help="run simulations (and their analysis)")
    p.add_argument("-p", "--plumed", nargs="+", required=True, help="PLUMED install label(s) in opt/plumed")
    p.add_argument("--set", help=f"name of the run set in {RUNS} (default: the PLUMED label)")
    add_selection(p)
    p.add_argument("-r", "--replicas", help="replicas to run, e.g. 0-4,7 (default: all in the tier)")
    p.add_argument("-j", "--jobs", type=int, default=os.cpu_count(), help="parallel jobs (default: all cores)")
    p.add_argument("--gromacs", default=DEFAULT_GROMACS, help=f"GROMACS install (default: {DEFAULT_GROMACS})")
    p.add_argument("--force", action="store_true", help="rerun also completed runs")
    p.add_argument("--dry-run", action="store_true", help="only print what would run")
    p.set_defaults(func=cmd_run)

    p = sub.add_parser("analyze", help="redo the analysis of completed runs")
    p.add_argument("-p", "--plumed", nargs="+", required=True, help="PLUMED install label(s)")
    p.add_argument("--set", help="name of the run set (default: the PLUMED label)")
    add_selection(p)
    p.add_argument("-r", "--replicas")
    p.add_argument("--gromacs", default=DEFAULT_GROMACS)
    p.set_defaults(func=cmd_analyze)

    p = sub.add_parser("compare", help="compare two run sets and write a report")
    p.add_argument("A", help="reference run set")
    p.add_argument("B", help="candidate run set")
    add_selection(p)
    p.add_argument("-o", "--out", help="output directory (default: reports/<A>_vs_<B>/<tier>)")
    p.add_argument("--strict", action="store_true", help="exit with error if any metric is marked")
    p.set_defaults(func=cmd_compare)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
