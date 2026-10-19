"""Command-line interface for IMEXscan."""

import argparse
import os
import shutil
import sys

from . import __version__
from . import dif as dif_mod
from . import pipeline


def _require_executable(executable, hint):
    if shutil.which(executable) is None:
        sys.exit(f"ERROR: '{executable}' not found. Activate the environment that "
                 f"provides it (e.g. conda activate imexscan) or {hint}.")


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="imexscan",
        description="IMEXscan: mining filamentous (Inovirus) prophages and their "
                    "pHVR regions from bacterial genomes")
    parser.add_argument("-v", "--version", action="version",
                        version=f"%(prog)s {__version__}")
    sub = parser.add_subparsers(dest="command", required=True, metavar="{run,dif,defense}")

    p_run = sub.add_parser("run", help="run the full mining pipeline on a directory of genomes")
    p_run.add_argument("-i", "--input", required=True, metavar="DIR",
                       help="directory with genome FASTA files (.fa/.fasta/.fna, optionally .gz)")
    p_run.add_argument("-o", "--output", required=True, metavar="DIR",
                       help="output directory")
    p_run.add_argument("-t", "--threads", type=int, default=4,
                       help="number of genomes processed in parallel (default: 4)")
    p_run.add_argument("--hmm", metavar="DIR", default=None,
                       help="marker HMM directory (default: bundled Vibrio inovirus HMMs)")
    p_run.add_argument("--min-hallmark", type=int, default=4,
                       help="minimum inovirus hallmark signals per prophage (default: 4)")
    p_run.add_argument("--evalue", default="1e-5", metavar="FLOAT",
                       help="hmmsearch E-value cutoff (default: 1e-5)")
    p_run.add_argument("--skip-dif", action="store_true",
                       help="skip the dif-like att motif scan (boundaries then "
                            "default to core +/- 5 kb, quality becomes 'incomplete')")
    p_run.add_argument("--prodigal", default="prodigal", metavar="PATH",
                       help="prodigal executable (default: prodigal from PATH)")
    p_run.add_argument("--hmmsearch", default="hmmsearch", metavar="PATH",
                       help="hmmsearch executable (default: hmmsearch from PATH)")
    p_run.set_defaults(func=_cmd_run)

    p_dif = sub.add_parser("dif", help="only scan a genome for dif-like att motifs")
    p_dif.add_argument("-i", "--input", required=True, metavar="FASTA")
    p_dif.add_argument("-o", "--output", default=None, metavar="CSV",
                       help="output CSV (default: <stem>-find_Vibrio_dif_motif-output.csv)")
    p_dif.set_defaults(func=_cmd_dif)

    p_df = sub.add_parser("defense",
                          help="annotate pHVR proteins with DefenseFinder "
                               "(requires defense-finder in PATH)")
    p_df.add_argument("-i", "--input", required=True, metavar="DIR",
                      help="output directory of `imexscan run`")
    p_df.add_argument("-o", "--output", required=True, metavar="DIR")
    p_df.add_argument("-t", "--threads", type=int, default=8,
                      help="threads per defense-finder run (default: 8)")
    p_df.add_argument("--chunk-size", type=int, default=5000,
                      help="proteins per chunk (default: 5000)")
    p_df.set_defaults(func=_cmd_defense)

    args = parser.parse_args(argv)
    args.func(args)


def _cmd_run(args):
    _require_executable(args.prodigal, "pass --prodigal PATH")
    _require_executable(args.hmmsearch, "pass --hmmsearch PATH")
    pipeline.run(input_dir=args.input, out_dir=args.output, threads=args.threads,
                 hmm_dir=args.hmm, min_hallmark=args.min_hallmark, evalue=args.evalue,
                 skip_dif=args.skip_dif, prodigal=args.prodigal, hmmsearch=args.hmmsearch)


def _cmd_dif(args):
    records = dif_mod.scan(args.input)
    out = args.output
    if out is None:
        base = os.path.basename(args.input)
        if base.endswith(".gz"):
            base = base[:-3]
        out = os.path.splitext(base)[0] + "-find_Vibrio_dif_motif-output.csv"
    dif_mod.write_csv(records, out)
    print(f"{len(records)} dif-like matches -> {out}")


def _cmd_defense(args):
    _require_executable("defense-finder",
                        "install it with conda install -c bioconda defense-finder")
    pipeline.run_defense(args.input, args.output, threads=args.threads,
                         chunk_size=args.chunk_size)
