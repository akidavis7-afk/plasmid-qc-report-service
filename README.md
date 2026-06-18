# Plasmid Mutation QC Report Service

A Python-based portfolio project for plasmid sequence comparison, mutation reporting, and research-use visualization.

This project demonstrates a custom workflow where researchers provide reference plasmid files and sample/mutated sequence files, and I generate clean output reports such as mutation tables, alignment previews, plasmid map figures, and summary PDFs.

## Service Overview

Many molecular biology labs need to compare expected plasmid sequences with actual mutated or sequencing-result sequences. This workflow helps convert those files into clear, reviewable outputs.
This service focuses on accurate mutation detection, reference-coordinate reporting, and clean summary deliverables. The main visualization output is a stable linear mutation map, with optional circular plasmid map output when suitable.

The service can provide:

* Reference vs sample sequence comparison
* Detection of mismatches, insertions, and deletions
* Mutation positions reported using original plasmid coordinates
* CSV mutation tables
* Alignment preview files
* Simple plasmid map figures in PNG/PDF
* Final Markdown/PDF summary reports

## Visualization Status

The current workflow uses a **linear mutation map** as the main stable visualization output. This map clearly shows mutation positions based on the original/reference plasmid coordinates.

A **circular plasmid map** can also be included as an optional experimental visualization when suitable. However, circular plasmid maps can become crowded when many GenBank features and mutation labels overlap, so the linear map is used as the primary reliable output for mutation review.

This is not intended to replace professional tools such as SnapGene, Benchling, or Addgene tools. The goal is to provide custom automation and reporting support for repetitive plasmid QC workflows.

## Example Workflow

```text
Client provides:
1. Reference GenBank file
2. Reference FASTA file
3. Sample or mutated FASTA file
4. Plasmid topology: circular or linear
5. Important features to show, if any

I provide:
1. mutation_report.csv
2. alignment.txt
3. summary_report.md
4. summary_report.pdf
5. linear_mutation_map.png
6. linear_mutation_map.pdf
7. optional circular_plasmid_map.png / .pdf when suitable
```

## Example Deliverables

After processing the provided plasmid files, the output package may include:

```text
outputs/
├── tables/
│   └── mutation_report.csv
│
├── reports/
│   ├── alignment.txt
│   ├── summary_report.md
│   └── summary_report.pdf
│
└── figures/
    ├── linear_mutation_map.png
    ├── linear_mutation_map.pdf
    ├── circular_plasmid_map.png      # optional / experimental
    └── circular_plasmid_map.pdf      # optional / experimental

```

## Example Mutation Report

Example output table:

| ID | Type      | Reference position | Reference sequence | Sample sequence | Length |
| -: | --------- | ------------------ | ------------------ | --------------- | -----: |
|  1 | mismatch  | 245                | A                  | G               |      1 |
|  2 | deletion  | 500-506            | ATCGGAA            | -               |      7 |
|  3 | insertion | after 1200         | -                  | TTGACA          |      6 |

Mutation positions are reported using the original/reference plasmid coordinate system.

## Input Files Needed

For a basic report, please provide:

1. **Reference plasmid file**

   * Preferred: GenBank file `.gb` or `.gbk`
   * Accepted: FASTA file `.fasta` or `.fa`

2. **Sample or mutated sequence file**

   * Preferred: FASTA file `.fasta` or `.fa`

3. **Topology information**

   * Circular or linear

4. **Important features to show**

   * Example: gene, promoter, ori, antibiotic marker, MCS, tag, restriction site

5. **Expected mutation information, if known**

   * Example: “We expected A245G”
   * Example: “There should be a small insertion near AmpR”

## What I Deliver

Depending on the request, I can provide:

- `mutation_report.csv`
- `alignment.txt`
- `summary_report.md`
- `summary_report.pdf`
- `linear_mutation_map.png`
- `linear_mutation_map.pdf`

Optional visualization:

- `circular_plasmid_map.png`
- `circular_plasmid_map.pdf`

The linear mutation map is the main stable visual output. The circular plasmid map is optional and may be included when the GenBank annotations are suitable for clear visualization.

For custom requests, I can also prepare:

* batch reports for multiple plasmids
* custom table formats
* lab-notebook-friendly summaries
* presentation-ready figures
* reusable private workflow scripts

## What This Portfolio Demonstrates

This project demonstrates my ability to:

* Work with FASTA and GenBank files
* Use Python for sequence analysis automation
* Compare reference and sample sequences
* Detect mismatches, insertions, and deletions
* Report sequence changes using reference coordinates
* Generate CSV, Markdown, PDF, and image outputs
* Build repeatable research-use reporting workflows

## Intended Use

This service is intended for:

* research-use plasmid QC support
* automation of repetitive sequence comparison tasks
* report generation
* figure generation
* portfolio demonstration
* non-clinical research workflows

## Not Intended For

This service is not intended for:

* clinical diagnosis
* medical decision-making
* patient data analysis
* final biological interpretation
* replacing professional molecular biology software
* large-scale NGS analysis
* complex genome assembly
* confidential data processing without prior agreement

## Limitations

* Mutation positions are reported based on the original/reference plasmid coordinate system.
* Insertions and deletions in repetitive regions may have ambiguous exact positions depending on alignment scoring.
* Large insertions, deletions, inversions, or rearrangements may require manual review.
* The linear mutation map is the primary stable visualization output.
* Circular plasmid map generation is included only as an optional/experimental visualization.
* Circular plasmid visualization can become crowded when many features overlap.
* This workflow currently assumes reference and sample sequences are in the same orientation.
* Final biological interpretation should be reviewed by qualified researchers.

## Data Privacy

Please do not send patient data, clinical data, or highly confidential unpublished data without prior agreement.

For the first sample report, non-confidential plasmid examples are preferred.

## Code and Copyright Notice

This repository is provided as a portfolio demonstration of my plasmid mutation QC reporting workflow.

No open-source license is granted.

You may view this repository for evaluation purposes only. You may not copy, reuse, modify, redistribute, or use the source code in your own projects without written permission.

If you are interested in using this workflow, please contact me for a custom report or private service arrangement.

## Contact

If you have a non-confidential plasmid example and would like a sample report, please contact me with:

* reference GenBank or FASTA file
* sample/mutated FASTA file
* circular or linear topology
* important features to show
* preferred output format: CSV, PNG, PDF, or all
