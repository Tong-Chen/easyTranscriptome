#!/usr/bin/env python3
"""Generate or execute a repository-backed RNA-seq Nextflow command."""

from __future__ import annotations

import argparse
import os
import shlex
import subprocess
import sys
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("salmon", "star", "both"), required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="Root for relative FASTQ paths; defaults to the metadata directory.",
    )
    parser.add_argument("--genome-fa", type=Path, required=True)
    parser.add_argument("--genome-gtf", type=Path, required=True)
    parser.add_argument("--genome-version", required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--group-col", default="conditions")
    parser.add_argument("--with-de", action="store_true")
    parser.add_argument("--with-stringtie", action="store_true")
    parser.add_argument("--salmon-index", type=Path)
    parser.add_argument("--star-index", type=Path)
    parser.add_argument("--pipeline-root", type=Path)
    parser.add_argument("--work-dir", type=Path, default=Path("."))
    parser.add_argument("--execute", action="store_true")
    parser.add_argument("--no-resume", action="store_true")
    return parser.parse_args()


def pipeline_root(explicit: Path | None) -> Path:
    if explicit:
        return explicit.resolve()
    env_root = os.environ.get("TRANSCRIPTOME_PIPELINE_ROOT")
    if env_root:
        return Path(env_root).expanduser().resolve()
    here = Path(__file__).resolve()
    candidates = [here.parents[4] / "pipeline", Path.cwd() / "pipeline"]
    found = next((path for path in candidates if (path / "nextflow").is_dir()), None)
    if found is None:
        raise ValueError("pipeline repository not found; pass --pipeline-root")
    return found.resolve()


def main() -> int:
    args = parse_args()
    try:
        root = pipeline_root(args.pipeline_root)
        entry = {
            "salmon": "RNAseq_salmon.nf",
            "star": "RNAseq_star.nf",
            "both": "RNAseq_salmon_star.nf",
        }[args.mode]
        entry_path = root / "nextflow" / entry
        if not entry_path.is_file():
            raise ValueError(f"Nextflow entry does not exist: {entry_path}")
        if args.with_stringtie and args.mode == "salmon":
            raise ValueError("StringTie requires --mode star or --mode both")
        if not args.metadata.is_file():
            raise ValueError(f"metadata does not exist: {args.metadata}")
        data_dir = (args.data_dir or args.metadata.parent).resolve()
        if not data_dir.is_dir():
            raise ValueError(f"data directory does not exist: {data_dir}")

        preflight_command = [
            sys.executable,
            str(Path(__file__).with_name("preflight.py")),
            "--metadata",
            str(args.metadata.resolve()),
            "--data-dir",
            str(data_dir),
            "--genome-fa",
            str(args.genome_fa.resolve()),
            "--genome-gtf",
            str(args.genome_gtf.resolve()),
            "--mode",
            args.mode,
            "--pipeline-root",
            str(root),
        ]
        if args.with_de:
            preflight_command.append("--with-de")
        preflight_result = subprocess.run(preflight_command, check=False)
        if preflight_result.returncode != 0:
            raise ValueError("metadata/input preflight failed; no workflow command generated")

        command = [
            "nextflow",
            "run",
            str(entry_path),
        ]
        if not args.no_resume:
            command.append("-resume")
        command.extend([
            "--metadata",
            str(args.metadata.resolve()),
            "--INPUT_CHECK_data_dir",
            str(data_dir),
            "--genome_fa",
            str(args.genome_fa.resolve()),
            "--genome_gtf",
            str(args.genome_gtf.resolve()),
            "--genome_version",
            args.genome_version,
            "--outdir",
            str(args.outdir.resolve()),
            "--group_col",
            args.group_col,
            "--skip_de_analysis",
            str(not args.with_de).lower(),
        ])
        if args.mode in {"star", "both"}:
            command.extend(
                ["--skip_stringtie_assembl", str(not args.with_stringtie).lower()]
            )
        if args.salmon_index:
            if args.mode == "star":
                raise ValueError("--salmon-index is incompatible with --mode star")
            command.extend(["--salmon_index", str(args.salmon_index.resolve())])
        if args.star_index:
            if args.mode == "salmon":
                raise ValueError("--star-index is incompatible with --mode salmon")
            command.extend(["--star_index", str(args.star_index.resolve())])
        helper_dir = root / "bin" / "bash"
        env = os.environ.copy()
        env["PATH"] = str(helper_dir) + os.pathsep + env.get("PATH", "")
        rendered = "PATH=" + shlex.quote(str(helper_dir)) + os.pathsep + '"$PATH" '
        rendered += shlex.join(command)
        print(rendered)
        if not args.execute:
            print("Dry run only. Re-run with --execute after preflight and command review.")
            return 0
        args.work_dir.mkdir(parents=True, exist_ok=True)
        completed = subprocess.run(command, cwd=args.work_dir, env=env, check=False)
        return completed.returncode
    except (OSError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
