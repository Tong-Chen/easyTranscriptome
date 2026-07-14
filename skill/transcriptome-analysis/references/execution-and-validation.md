# Execution and validation

## Environment profiles

### Local computer

Use Salmon for routine work when CPU, RAM, or storage is limited. Build indexes once and reuse them. Run Nextflow from a dedicated project directory so `.nextflow.log`, `work/`, and result paths are isolated.

The current repository modules reference `/anaconda3/envs/transcriptome` and, for Salmon index construction, `/anaconda3/envs/agat`. If those paths do not exist, either adapt a copy of the modules to the local Conda installation or run on the configured server. Do not create system-level paths automatically.

### Server

Confirm scheduler policy, scratch space, quota, CPU/memory limits, Conda paths, and retention policy. The checked-in Nextflow config defaults to local execution; configure the server executor separately when required. Use `-resume` after a recoverable failure and preserve the original command and trace.

## Commands

Generate a Salmon command with differential analysis:

```bash
python3 scripts/run_transcriptome.py --mode salmon \
  --metadata metadata.tsv --data-dir raw_data_root \
  --genome-fa genome.fa --genome-gtf genes.gtf \
  --genome-version build_name --outdir result --group-col conditions --with-de
```

Generate a STAR command with StringTie but without differential analysis:

```bash
python3 scripts/run_transcriptome.py --mode star \
  --metadata metadata.tsv --data-dir raw_data_root \
  --genome-fa genome.fa --genome-gtf genes.gtf \
  --genome-version build_name --outdir result --with-stringtie
```

Reuse an index with `--salmon-index PATH` or `--star-index PATH`. Omit `--data-dir` only when relative FASTQ paths are based at the metadata directory or all paths are absolute. Generate the command without `--execute`, inspect it, and then explicitly execute after preflight.

## Expected result surfaces

- `result/pipeline_info/`: execution report, timeline, trace, and DAG.
- `result/salmon_index/` or `result/star_index/`: reusable indexes when built.
- `result/salmon_quant/`: per-sample Salmon quantification and summary table.
- `result/star_align/`: BAM/CSI, STAR logs/counts, and merged count matrix.
- `result/salmon_de/` or `result/readscount_de/`: normalized expression, PCA, correlation, DE tables, heatmaps, volcano/rank plots, and parameters.
- `result/stringtie/`: per-sample and merged transcript annotations.

Exact directories depend on the module `publishDir` declarations; verify against the execution trace rather than assuming success.

## Completion checklist

1. Nextflow returns zero and no process is failed or silently ignored.
2. Trace contains the expected number of sample processes.
3. Every metadata sample has exactly one nonempty `quant.sf` or expected STAR BAM/count file.
4. FASTQ and mapping/quantification QC are reviewed, not merely generated.
5. Count matrices contain the expected sample columns and identifier type.
6. PCA and correlation agree with the recorded design or discrepancies are explained.
7. DE comparison direction, design formula, covariates, log2FC threshold, and adjusted-P threshold are recorded.
8. The handoff distinguishes completed computation, statistical evidence, and biological interpretation.
