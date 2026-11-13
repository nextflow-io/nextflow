---
name: metadata-ontology
description: Map free-text sequencing-sample metadata to controlled-vocabulary ontology terms. Use this skill whenever labelling a sample's organism, tissue/cell type, disease, or assay from its metadata.
---
# Sample metadata → ontology labelling

You assign controlled-vocabulary labels to a sequencing sample from its (often
messy, free-text) metadata. Follow these rules exactly.

## Output convention

Every facet is either a single term written as `label [ONTOLOGY:ID]`, or the
literal string `unknown` when the metadata does not evidence it. Never invent an
ID: only use IDs from the tables below (or a `NCBITaxon:` id when the record
gives an explicit `tax_id`). If the right term is not in a table, use the
closest listed parent term, or `unknown`.

## Allowed terms (illustrative curated subset — replace with your own vocabulary)

**Organism** (from `scientific_name` / `tax_id`):
| term | id |
|---|---|
| Homo sapiens | NCBITaxon:9606 |
| Mus musculus | NCBITaxon:10090 |
| Drosophila melanogaster | NCBITaxon:7227 |
| Saccharomyces cerevisiae | NCBITaxon:4932 |

**Tissue / cell type** (UBERON tissues, CL cell types):
| term | id |
|---|---|
| B lymphocyte | CL:0000236 |
| epithelial cell | CL:0000066 |
| lung | UBERON:0002048 |
| blood | UBERON:0000178 |
| liver | UBERON:0002107 |
| breast | UBERON:0000310 |

**Disease** (MONDO):
| term | id |
|---|---|
| lung adenocarcinoma | MONDO:0005097 |
| none | — |

Use `none [—]` only when the sample is clearly a normal/healthy or engineered
non-disease sample; otherwise `unknown`.

**Assay** (functional-genomics assay term):
| term |
|---|
| RNA-seq |
| ATAC-seq |
| ChIP-seq |
| WGS |
| WES |

Write the assay as the bare term (no id).

## Rules

- **Read the free-text title, do not trust structured fields blindly.** The
  `library_strategy` field is frequently `OTHER` or wrong; if the `sample_title`
  says `ATAC`/`ATACseq`, the assay is `ATAC-seq` regardless of `library_strategy`.
- **Recognize common cell lines** and label their cell type: `GM12878` is a
  lymphoblastoid line → `B lymphocyte [CL:0000236]`; `A549` → lung epithelial;
  `S2`/`S2-DRSC` is a *Drosophila* embryonic line.
- **Set confidence to the evidence.** A facet named explicitly in the title →
  high; inferred from a cell-line name → medium; not evidenced at all → output
  `unknown` and lower the overall confidence.
- **Never fabricate.** A bare replicate label like `N61311_untreated` or
  `Set7KD_rep1` evidences no tissue or disease — return `unknown`, not a guess.
- In the rationale, quote the exact phrase(s) you keyed on, and name any facet
  you left `unknown`.
