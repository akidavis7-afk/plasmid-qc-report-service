# SYNPUC19CV mutation QC report

Generated: 2026-06-17 19:45

## Overview

This report compares a reference plasmid sequence with a mutated or sample sequence and summarizes detected sequence differences.

This report is intended for research-use automation and visualization. Final biological interpretation should be reviewed by the researcher.

## Input files

| Item | Path |
|---|---|
| Reference GenBank | `data/reference/reference_puc19.gb` |
| Mutation report CSV | `outputs/tables/mutation_report.csv` |
| Alignment preview | `outputs/reports/alignment.txt` |
| Plasmid map image | `outputs/figures/puc19_mutation_map.png` |

## Plasmid information

| Field | Value |
|---|---|
| Plasmid name | SYNPUC19CV |
| Length | 2,686 bp |
| Topology | circular |
| GenBank features | 1 |

## Mutation summary

| Type | Count |
|---|---:|
| Total mutations | 3 |
| Mismatches | 1 |
| Insertions | 1 |
| Deletions | 1 |
| Other | 0 |

## Detected mutations

| ID | Type | Reference position | Reference sequence | Sample sequence | Length | Note |
|---:|---|---|---|---|---:|---|
| 1 | mismatch | 245 | G | A | 1 | G changed to A at reference position 245 |
| 2 | deletion | 500-506 | TCACAAT | - | 7 | 7 bp deleted from reference position 500-506 |
| 3 | insertion | after 1200 | - | TTGACA | 6 | 6 bp inserted after 1200 |


## Selected GenBank features

No selected GenBank features found.


## Output files

| File | Description |
|---|---|
| `mutation_report.csv` | Table of detected sequence changes |
| `alignment.txt` | Gapped reference/sample alignment preview |
| `summary_report.md` | Markdown summary report |
| `summary_report.pdf` | PDF summary report |
| plasmid map PNG/PDF | Visual map generated separately by Script 4 |

## Limitations

- Mutation positions are reported using the original/reference plasmid coordinate system.
- Insertions and deletions in repetitive regions may have ambiguous exact positions depending on alignment scoring.
- This script does not provide clinical or diagnostic interpretation.
- Large rearrangements may require manual review or a more advanced alignment workflow.
