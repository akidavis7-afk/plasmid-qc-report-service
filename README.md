# Plasmid QC Report Service

A portfolio demonstration of a Python-based plasmid mutation QC reporting workflow.

Researchers provide plasmid sequence files, and I provide clean output reports such as mutation tables, alignment previews, summary PDFs, and simple mutation map figures.

This service focuses on **custom reporting and repetitive workflow automation**. It does not replace professional plasmid software such as SnapGene, Benchling, Addgene tools, or sequencing providers.

---

## What This Service Does

This workflow compares a reference plasmid sequence with a sample or mutated sequence and reports:

* mismatches
* insertions
* deletions
* mutation positions based on the original reference plasmid coordinates
* alignment preview
* simple mutation map figure
* summary report in Markdown/PDF format

The main stable visualization is a **linear mutation map**. A circular plasmid map may be included as an optional experimental visualization when suitable.

---

## What the Client Provides

For a basic report, please provide:

1. Reference plasmid file

   * Preferred: GenBank `.gb` / `.gbk`
   * Accepted: FASTA `.fasta` / `.fa`

2. Sample or mutated sequence file

   * Preferred: FASTA `.fasta` / `.fa`

3. Plasmid topology

   * circular or linear

4. Important features to show, if any

   * gene
   * promoter
   * ori
   * antibiotic marker
   * MCS
   * tag
   * restriction site

5. Expected mutation information, if known

   * Example: “We expected A245G”
   * Example: “There should be an insertion near AmpR”

For a first sample report, non-confidential plasmid examples are preferred.

---

## What I Deliver

Depending on the request, I can provide:

* `mutation_report.csv`
* `alignment_preview_sample.txt`
* `summary_report.md`
* `summary_report.pdf`
* `linear_mutation_map.png`
* `linear_mutation_map.pdf`

Optional output:

* `circular_plasmid_map.png`
* `circular_plasmid_map.pdf`

The linear mutation map is the main stable visual output. Circular plasmid maps can become crowded when many GenBank features and mutation labels overlap, so circular output is optional.

---

## Example Deliverables

```text
examples/
├── mutation_report_sample.csv
├── alignment_preview_sample.txt
├── summary_report_sample.md
├── summary_report_sample.pdf
└── figures/
    ├── linear_mutation_map_sample.png
    └── linear_original_map_sample.png
```

Example mutation report:

| ID | Type      | Reference position | Reference sequence | Sample sequence | Length |
| -: | --------- | ------------------ | ------------------ | --------------- | -----: |
|  1 | mismatch  | 245                | A                  | G               |      1 |
|  2 | deletion  | 500-506            | ATCGGAA            | -               |      7 |
|  3 | insertion | after 1200         | -                  | TTGACA          |      6 |

Mutation positions are reported using the original/reference plasmid coordinate system.

---

## Intended Use

This service is intended for:

* research-use plasmid QC support
* repetitive sequence comparison reporting
* mutation table generation
* summary PDF generation
* simple figure generation
* non-clinical research workflows

---

## Limitations

* This service is for research-use automation and visualization only.
* It does not provide clinical, diagnostic, or final biological interpretation.
* Final biological interpretation should be reviewed by qualified researchers.
* Mutation positions are reported based on the original/reference plasmid coordinate system.
* Insertions and deletions in repetitive regions may have ambiguous exact positions depending on alignment scoring.
* Large insertions, deletions, inversions, or rearrangements may require manual review.
* Circular plasmid map generation is optional and may not be suitable for highly crowded GenBank annotations.
* This workflow assumes the reference and sample sequences are in the same orientation.

---

## Data Privacy

Please do not send patient data, clinical data, or highly confidential unpublished data without prior agreement.

For a first sample report, please use non-confidential plasmid files.

---

## Code Availability

This repository is a public service portfolio. It shows example deliverables and service information.

The full working Python workflow is not included in this public repository.

No open-source license is granted. The workflow, templates, and source code are not available for copying, reuse, modification, or redistribution without written permission.

---

## Contact

To request a sample report, please send:

* reference GenBank or FASTA file
* sample/mutated FASTA file
* circular or linear topology
* important features to show
* preferred output format: CSV, PNG, PDF, or all
