# Input contract

## Required user-supplied metadata

Require metadata as an input to every analysis. Do not construct it by scanning a FASTQ directory: file naming conventions do not reliably encode sample identity, mate assignment, biological condition, batch, or comparison direction.

Use UTF-8, tab separation, one row per biological sample, and this exact minimum header:

```text
sample	conditions	single_end	fastq_1	fastq_2
WT_1	WT	false	seq/WT_1_R1.fastq.gz	seq/WT_1_R2.fastq.gz
WT_2	WT	false	seq/WT_2_R1.fastq.gz	seq/WT_2_R2.fastq.gz
KO_1	KO	false	seq/KO_1_R1.fastq.gz	seq/KO_1_R2.fastq.gz
KO_2	KO	false	seq/KO_2_R1.fastq.gz	seq/KO_2_R2.fastq.gz
```

- `sample`: unique, stable, shell-safe sample identifier.
- `conditions`: biological group used by DESeq2; do not encode batch here.
- `single_end`: literal `true` or `false`.
- `fastq_1`: read 1 or the single-end FASTQ.
- `fastq_2`: read 2 for paired data; empty only for single-end data.

Require every declared FASTQ to exist and contain a readable FASTQ record before running. Resolve absolute paths directly. Resolve every relative path against one explicit data root supplied as `--data-dir`; when omitted, use the metadata file's directory. Pass the same root to Nextflow as `INPUT_CHECK_data_dir` so validation and execution cannot disagree.

Reject metadata when a paired sample lacks either mate, both mate columns point to the same file, a single-end row supplies `fastq_2`, or the same FASTQ is assigned to more than one sample/read slot.

## Design and comparisons

Keep contrast direction explicit, for example `KO` versus `WT`, because log2 fold-change signs depend on numerator/reference ordering. Put batch, sex, donor, or paired-subject effects in separate metadata columns and pass the appropriate covariate option when the downstream wrapper supports it.

## Design checks

- Require independent biological replicates. Technical FASTQ lanes from the same library are not biological replicates; concatenate or model them correctly before sample-level DE.
- Prefer at least three biological replicates per group. Flag one replicate as invalid for dispersion-based differential inference and two as fragile.
- Do not mix genome builds, annotation releases, gene identifier namespaces, or incompatible chromosome naming conventions.
- Record strandedness from the library protocol; do not guess it from paired-end status.
- Record exclusions with reasons and retain the original sample table.
- Avoid spaces, tabs, slashes, quotes, and shell metacharacters in sample IDs.
