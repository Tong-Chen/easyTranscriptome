---
name: transcriptome-analysis
description: Plan, prepare, run, and validate bulk RNA-seq transcriptome analyses from a required user-supplied metadata table using the easyTranscriptome local/server teaching workflows and the reusable Nextflow pipelines under pipeline/nextflow. Use when Codex needs to validate that metadata resolves every declared FASTQ, choose Salmon versus STAR, generate a reproducible command, perform expression quantification or DESeq2 differential analysis, assemble transcripts with StringTie, troubleshoot a run, or summarize analysis outputs on a local computer or server.
---

# Transcriptome analysis

Use the repository's existing workflows as the source of truth. Prefer the DSL2 Nextflow entry points for reproducible batch execution; use the local R Markdown files for downstream interpretation and the server shell workflow for method details and recovery guidance.

## Resolve repository paths

Find the skill directory first, then resolve these sibling repositories without assuming the current working directory:

- `easyTranscriptome`: three directories above this file.
- Nextflow root: `$TRANSCRIPTOME_PIPELINE_ROOT/nextflow` when the environment variable exists; otherwise try the sibling repository `pipeline/nextflow` under the common parent.
- Pipeline helper commands: `$TRANSCRIPTOME_PIPELINE_ROOT/bin/bash`.

Do not edit raw FASTQ, reference FASTA/GTF, or an existing result directory in place. Create a new output directory or use Nextflow `-resume` for an interrupted run.

## Choose the analysis route

Read [references/workflow-map.md](references/workflow-map.md) when choosing stages or tracing an implementation back to the original files.

- Choose `salmon` for fast transcript/gene quantification and routine differential expression when BAM files, novel transcripts, and splice events are unnecessary.
- Choose `star` when coordinate-sorted BAM, genome-browser inspection, gene counts, StringTie assembly, or alternative-splicing analysis is required.
- Choose `both` only when the user needs alignment artifacts plus Salmon quantification or a cross-method comparison; explain the additional compute and storage.
- Treat de novo assembly, novel lncRNA, ceRNA, rMATS, WGCNA, GSEA, and enrichment as explicit extensions. Do not silently add them to a standard expression workflow.

## Prepare inputs

Read [references/input-contract.md](references/input-contract.md) before validating metadata.

1. Require the user to provide metadata. Do not generate it from FASTQ filenames or infer sample identity, pairing, biological groups, covariates, or comparisons from directory contents.
2. Inventory reference FASTA/GTF, library layout, strandedness, organism, genome build, and the data root used by relative FASTQ paths.
3. Preserve all user-supplied sample IDs and design fields.
4. Run preflight before constructing an execution command. Resolve absolute FASTQ paths directly; resolve relative paths only against `--data-dir`, which defaults to the metadata directory:

   ```bash
   python3 scripts/preflight.py --metadata metadata.tsv --data-dir raw_data_root \
     --genome-fa genome.fa --genome-gtf genes.gtf --mode salmon --with-de
   ```

Stop if metadata is absent or preflight reports any unresolved, unreadable, invalid, duplicated, or inconsistent FASTQ. Surface warnings about low replication, `NA` conditions, missing tools, hard-coded environment paths, or mismatched references before computing.

## Build and run the workflow

Read [references/execution-and-validation.md](references/execution-and-validation.md) for environment constraints, command examples, output locations, and completion checks.

Generate a command first:

```bash
python3 scripts/run_transcriptome.py \
  --mode salmon \
  --metadata metadata.tsv \
  --data-dir raw_data_root \
  --genome-fa genome.fa \
  --genome-gtf genes.gtf \
  --genome-version GRCh38 \
  --outdir result \
  --group-col conditions \
  --with-de
```

Inspect the printed command and preflight report. The runner repeats preflight and refuses to generate or execute a command when metadata cannot resolve its FASTQ inputs. Add `--execute` only when the user requested execution and the selected machine has adequate compute, storage, software, and the expected Conda paths. Keep the generated command in the handoff even when execution is not authorized.

For differential analysis, require a real condition column with at least two groups and preferably at least three biological replicates per group. Capture covariates and comparison pairs explicitly. Do not claim causal biological conclusions from the pipeline output alone.

## Validate and report

Validate in this order:

1. Confirm Nextflow exit status and inspect `pipeline_info` report, timeline, trace, and DAG.
2. Confirm every metadata sample has its expected quantification or alignment output.
3. Review FastQC/MultiQC and Salmon mapping rate or STAR uniquely mapped/multimapped/unmapped metrics before interpreting counts.
4. For DESeq2, inspect library size, sample correlation, PCA, normalized-expression distribution, dispersion/outliers, comparison direction, adjusted P value, and fold-change thresholds.
5. Report exact input build, sample exclusions, parameters, software route, warnings, failed/skipped stages, output paths, and what remains biologically unverified.

Never describe a run as complete solely because output files exist. Require integrity checks and one-to-one coverage of all metadata samples.
