# Workflow map

## Route selection

| Need | Route | Primary implementation | Main outputs |
| --- | --- | --- | --- |
| Fast gene/transcript abundance | Salmon | `pipeline/nextflow/RNAseq_salmon.nf` | per-sample `quant.sf`, Salmon output table |
| Genome alignment and gene count | STAR | `pipeline/nextflow/RNAseq_star.nf` | sorted BAM/CSI, `ReadsPerGene.out.tab`, count matrix |
| Both quantification strategies | Both | `pipeline/nextflow/RNAseq_salmon_star.nf` | Salmon outputs plus STAR outputs |
| Differential expression | Add `--with-de` | `Salmon_quant_DE.nf` or `Readscount_DE.nf` | normalized matrices, PCA/correlation, DE tables and plots |
| Novel transcript assembly | STAR plus `--with-stringtie` | `Stringtie_assemble.nf` | per-sample GTF and merged/annotated GTF |
| Alternative splicing | STAR extension | `server/RNAseq_pipeline.sh`, sections 4–5 | rMATS event tables and sashimi plots |
| De novo transcriptome | Explicit extension | `server/RNAseq_pipeline.sh`, section 6 | Trinity assembly and derived annotations |
| Novel lncRNA/ceRNA | Explicit extension | `server/RNAseq_pipeline.sh`, section 7 | coding-potential, target, correlation/network files |
| Local DE interpretation | Local R | `local_computer/13_salmon.DESeq2.simpler.Rmd`, `22_STAR.DESeq2.Batch.Rmd` | statistical tables and diagnostic plots |
| GO/GSEA/WGCNA | Downstream extension | `local_computer/14_GO_enrichment.Rmd`, `15_GO_visual.Rmd`, `16_GSEA/`, `17_WGCNA_simplest.Rmd` | enrichment and network results |

## Existing Nextflow module graph

- Input: `Input_check.nf` parses a tab-separated samplesheet into sample metadata and FASTQ channels.
- Salmon: `Salmon_build_index.nf` -> `Salmon_quant.nf` -> optional `Salmon_quant_DE.nf`.
- STAR: `STAR_build_index.nf` -> `STAR_align.nf` -> optional `Readscount_DE.nf`.
- Assembly: STAR coordinate-sorted BAM -> `Stringtie_assemble.nf`.
- Configuration: `nextflow.config` includes `RNAseq.config` and defines execution resources, reports, trace, timeline, and DAG.

## Boundaries

- The server shell file is a literate training workflow, not a parameterized production entry point. Reuse its commands selectively; do not execute the whole file blindly.
- The Windows/local shell files contain machine-specific paths and teaching examples. Treat them as method references.
- The current Nextflow modules call helper scripts from `pipeline/bin/bash` and contain absolute Conda environment paths. Preflight the target machine before execution.
- Reference FASTA, GTF, indexes, transcript-to-gene mapping, and gene identifier namespace must describe the same genome build.
