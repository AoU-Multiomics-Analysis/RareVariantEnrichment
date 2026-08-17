# Multi-Omics Phenotype ID Normalization Design

## Goal

Extend the LoF/PC rare-variant enrichment analysis to accept expression, protein-expression, and splicing outlier BED files while preserving the existing gene-level enrichment model.

## Observed phenotype ID formats

The Susie merged files in `~/Desktop/susie_files` use the following molecular phenotype identifiers:

- Expression: `ENSG00000000419.14`
- Proteomics: `A0JNW5_ENSG00000111647.13`
- Splicing: `chr20:50941209:50942031:clu_63027_-:ENSG00000000419.14`

All three formats contain an Ensembl gene identifier with an optional numeric version suffix. The existing carrier table and protein-coding gene table use unversioned Ensembl gene identifiers.

## Design

### Phenotype ID canonicalization

Add a phenotype-specific canonicalization function to the LoF/PC implementation. It will:

1. Strip surrounding whitespace.
2. Extract the Ensembl gene token matching `ENSG` followed by digits, with an optional numeric version suffix.
3. Remove the version suffix.
4. Reject empty or unsupported/ambiguous phenotype IDs with a line-numbered error.

The existing generic version-stripping behavior used for carrier and GTF gene IDs remains unchanged. No phenotype-type input is added to the WDL because the canonicalization is format-agnostic.

### Gene-level aggregation

The phenotype BED will be read into one gene-level phenotype vector per canonical Ensembl gene. If multiple input rows map to the same gene, each sample’s value will be the minimum finite z-score across those rows. Thus, a sample is represented by its most extreme negative splice-junction outlier for that gene.

Missing values are ignored during the minimum operation. A gene-sample value remains missing only when all contributing phenotype rows are missing. Expression and protein data with one row per gene are unchanged by this operation.

The downstream residualization, outlier thresholding, LoF carrier matching, and Fisher exact tests remain gene-level and otherwise unchanged.

### Quality control

The analysis QC JSON will distinguish input feature rows from normalized gene counts and report duplicate feature rows collapsed during aggregation. It will retain the existing counts for protein-coding genes and shared samples, adding explicit counts where needed rather than silently changing the meaning of existing fields.

### Workflow and outputs

The WDL public inputs and task boundaries remain unchanged. The same `lof-pc-enrichment` command will accept all three phenotype classes. Existing output files remain the interface:

- enrichment results TSV;
- summary JSON;
- gene/PC QC TSV;
- analysis QC JSON.

### Testing

Tests will cover:

- canonicalization of expression, proteomics, and splicing IDs using the observed Susie formats;
- rejection of unsupported or ambiguous IDs;
- per-sample minimum-z aggregation across multiple splice rows;
- missing-value behavior when some or all duplicate rows are missing;
- end-to-end enrichment with the existing fixture inputs;
- unchanged WDL input and output contracts.
