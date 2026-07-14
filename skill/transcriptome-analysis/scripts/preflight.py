#!/usr/bin/env python3
"""Validate RNA-seq inputs and runtime assumptions before expensive computation."""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import shutil
import sys
from collections import Counter
from pathlib import Path


REQUIRED_COLUMNS = {"sample", "conditions", "single_end", "fastq_1", "fastq_2"}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--genome-fa", type=Path, required=True)
    parser.add_argument("--genome-gtf", type=Path, required=True)
    parser.add_argument("--mode", choices=("salmon", "star", "both"), required=True)
    parser.add_argument("--with-de", action="store_true")
    parser.add_argument(
        "--data-dir",
        type=Path,
        help="Root for relative FASTQ paths; defaults to the metadata directory.",
    )
    parser.add_argument("--pipeline-root", type=Path)
    parser.add_argument("--json", action="store_true", help="Emit a JSON report.")
    return parser.parse_args()


def readable_file(path: Path, label: str, errors: list[str]) -> None:
    if not path.is_file():
        errors.append(f"{label} does not exist: {path}")
    elif path.stat().st_size == 0:
        errors.append(f"{label} is empty: {path}")


def fastq_has_record(path: Path) -> bool:
    opener = gzip.open if path.name.lower().endswith(".gz") else open
    try:
        with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
            lines = [handle.readline().rstrip("\n\r") for _ in range(4)]
        return len(lines) == 4 and lines[0].startswith("@") and lines[2].startswith("+")
    except (OSError, EOFError):
        return False


def resolve_read(raw: str, data_dir: Path) -> Path:
    path = Path(raw)
    if path.is_absolute():
        return path
    return data_dir / path


def locate_pipeline_root(explicit: Path | None) -> Path | None:
    if explicit:
        return explicit
    here = Path(__file__).resolve()
    candidates = [here.parents[4] / "pipeline", Path.cwd() / "pipeline"]
    return next((path for path in candidates if (path / "nextflow").is_dir()), None)


def main() -> int:
    args = parse_args()
    errors: list[str] = []
    warnings: list[str] = []
    data_dir = (args.data_dir or args.metadata.parent).resolve()
    readable_file(args.metadata, "metadata", errors)
    readable_file(args.genome_fa, "genome FASTA", errors)
    readable_file(args.genome_gtf, "genome GTF", errors)
    rows: list[dict[str, str]] = []

    if not errors:
        try:
            with args.metadata.open(newline="", encoding="utf-8-sig") as handle:
                reader = csv.DictReader(handle, delimiter="\t")
                missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
                if missing:
                    errors.append("metadata is missing columns: " + ", ".join(sorted(missing)))
                else:
                    rows = list(reader)
        except OSError as exc:
            errors.append(str(exc))

    samples: set[str] = set()
    seen_reads: dict[Path, tuple[str, str]] = {}
    conditions: Counter[str] = Counter()
    for line_number, row in enumerate(rows, 2):
        sample = (row.get("sample") or "").strip()
        condition = (row.get("conditions") or "").strip()
        single_raw = (row.get("single_end") or "").strip().lower()
        if not sample:
            errors.append(f"metadata:{line_number}: empty sample")
        elif sample in samples:
            errors.append(f"metadata:{line_number}: duplicate sample {sample}")
        samples.add(sample)
        if condition:
            conditions[condition] += 1
        if single_raw not in {"true", "false"}:
            errors.append(f"metadata:{line_number}: single_end must be true or false")
        if single_raw == "true" and (row.get("fastq_2") or "").strip():
            errors.append(f"metadata:{line_number}: fastq_2 must be empty for single-end data")
        read_fields = ["fastq_1"] if single_raw == "true" else ["fastq_1", "fastq_2"]
        row_paths: list[Path] = []
        for field in read_fields:
            raw = (row.get(field) or "").strip()
            if not raw:
                errors.append(f"metadata:{line_number}: {field} is required")
                continue
            path = resolve_read(raw, data_dir).resolve()
            row_paths.append(path)
            if path in seen_reads:
                previous_sample, previous_field = seen_reads[path]
                errors.append(
                    f"metadata:{line_number}: {field} reuses FASTQ already assigned to "
                    f"{previous_sample}.{previous_field}: {path}"
                )
            else:
                seen_reads[path] = (sample, field)
            if not path.is_file():
                errors.append(f"metadata:{line_number}: FASTQ does not exist: {path}")
            elif not fastq_has_record(path):
                errors.append(f"metadata:{line_number}: invalid or unreadable FASTQ: {path}")
        if len(row_paths) == 2 and row_paths[0] == row_paths[1]:
            errors.append(f"metadata:{line_number}: fastq_1 and fastq_2 resolve to the same file")

    if not rows and not errors:
        errors.append("metadata has no sample rows")
    if args.with_de:
        usable = {key: value for key, value in conditions.items() if key and key.upper() != "NA"}
        if len(usable) < 2:
            errors.append("DE analysis requires at least two non-NA conditions")
        for condition, count in sorted(usable.items()):
            if count < 2:
                warnings.append(f"condition {condition} has only {count} biological replicate")
            elif count < 3:
                warnings.append(f"condition {condition} has fewer than 3 biological replicates")

    pipeline_root = locate_pipeline_root(args.pipeline_root)
    if pipeline_root is None:
        errors.append("pipeline repository not found; pass --pipeline-root or set it in the runner")
    else:
        entry_names = {
            "salmon": ["RNAseq_salmon.nf"],
            "star": ["RNAseq_star.nf"],
            "both": ["RNAseq_salmon_star.nf"],
        }[args.mode]
        for name in entry_names:
            readable_file(pipeline_root / "nextflow" / name, "Nextflow entry", errors)
        if not (pipeline_root / "bin" / "bash").is_dir():
            warnings.append(f"pipeline helper directory not found: {pipeline_root / 'bin' / 'bash'}")

    required_tools = ["nextflow"]
    for tool in required_tools:
        if shutil.which(tool) is None:
            warnings.append(f"executable not on PATH: {tool}")
    if Path("/anaconda3/envs/transcriptome").exists() is False:
        warnings.append(
            "existing Nextflow modules reference /anaconda3/envs/transcriptome; adapt the modules or run on the configured server"
        )
    if args.mode in {"salmon", "both"} and not Path("/anaconda3/envs/agat").exists():
        warnings.append(
            "Salmon index construction references /anaconda3/envs/agat; provide a prebuilt index or adapt the module"
        )

    report = {
        "ok": not errors,
        "mode": args.mode,
        "samples": len(rows),
        "conditions": dict(sorted(conditions.items())),
        "data_dir": str(data_dir),
        "pipeline_root": str(pipeline_root) if pipeline_root else None,
        "errors": errors,
        "warnings": warnings,
    }
    if args.json:
        print(json.dumps(report, ensure_ascii=False, indent=2))
    else:
        print(f"Preflight: {'PASS' if report['ok'] else 'FAIL'}")
        print(f"Samples: {report['samples']}; conditions: {report['conditions']}")
        for warning in warnings:
            print(f"WARNING: {warning}")
        for error in errors:
            print(f"ERROR: {error}")
    return 0 if report["ok"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
