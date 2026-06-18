"""
generate_summary.py

Script 5 for plasmid-mutation-qc.

This script generates a clean Markdown and PDF summary report from:

1. Reference GenBank file
2. mutation_report.csv
3. Optional alignment.txt
4. Optional plasmid map PNG

Outputs:
- summary_report.md
- summary_report.pdf

This report is intended for research-use automation and visualization.
Final biological interpretation should be reviewed by the researcher.
"""

import argparse
import os
from pathlib import Path
from datetime import datetime

import pandas as pd
from Bio import SeqIO


# ------------------------------------------------------------
# Basic readers
# ------------------------------------------------------------

def read_genbank_record(genbank_path):
    if not os.path.exists(genbank_path):
        raise FileNotFoundError(f"GenBank file not found: {genbank_path}")

    records = list(SeqIO.parse(genbank_path, "genbank"))

    if not records:
        raise ValueError(f"No records found in GenBank file: {genbank_path}")

    return records[0]


def read_mutation_report(csv_path):
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"Mutation report not found: {csv_path}")

    df = pd.read_csv(csv_path)

    required_columns = [
        "mutation_id",
        "event_type",
        "reference_position",
        "reference_start",
        "reference_end",
        "reference_sequence",
        "sample_sequence",
        "length",
        "note",
    ]

    missing = [col for col in required_columns if col not in df.columns]

    if missing:
        raise ValueError(f"Mutation report is missing columns: {', '.join(missing)}")

    return df


def read_alignment_preview(alignment_path, max_chars=1200):
    """
    Read a short preview of alignment.txt.
    Avoid loading too much into the PDF.
    """
    if not alignment_path or not os.path.exists(alignment_path):
        return None

    with open(alignment_path, "r", encoding="utf-8") as file:
        content = file.read()

    if len(content) > max_chars:
        return content[:max_chars] + "\n\n... alignment preview truncated ..."

    return content


# ------------------------------------------------------------
# Data summarizers
# ------------------------------------------------------------

def get_plasmid_name(record):
    if record.name and record.name != "<unknown name>":
        return record.name

    if record.id and record.id != "<unknown id>":
        return record.id

    if record.description:
        return record.description.split(",")[0]

    return "Unknown plasmid"


def get_top_features(record, max_features=15):
    """
    Extract important GenBank features for summary report.
    """
    important_types = {
        "gene",
        "CDS",
        "promoter",
        "rep_origin",
        "ori",
        "misc_feature",
        "primer_bind",
        "terminator",
    }

    rows = []

    for feature in record.features:
        if feature.type not in important_types:
            continue

        try:
            start = int(feature.location.start) + 1
            end = int(feature.location.end)

            label = None

            for key in ["gene", "label", "product", "note"]:
                if key in feature.qualifiers and feature.qualifiers[key]:
                    label = str(feature.qualifiers[key][0])
                    break

            if not label:
                label = feature.type

            label = label.replace("\n", " ").strip()

            if len(label) > 40:
                label = label[:37] + "..."

            rows.append(
                {
                    "type": feature.type,
                    "label": label,
                    "start": start,
                    "end": end,
                    "strand": feature.strand if feature.strand else ".",
                }
            )

            if len(rows) >= max_features:
                break

        except Exception:
            continue

    return rows


def summarize_mutations(mutations_df):
    if mutations_df.empty:
        return {
            "total": 0,
            "mismatch": 0,
            "insertion": 0,
            "deletion": 0,
            "other": 0,
        }

    counts = mutations_df["event_type"].str.lower().value_counts().to_dict()

    return {
        "total": len(mutations_df),
        "mismatch": counts.get("mismatch", 0),
        "insertion": counts.get("insertion", 0),
        "deletion": counts.get("deletion", 0),
        "other": len(mutations_df)
        - counts.get("mismatch", 0)
        - counts.get("insertion", 0)
        - counts.get("deletion", 0),
    }


def shorten_text(value, max_length=60):
    value = str(value)

    if len(value) <= max_length:
        return value

    return value[: max_length - 3] + "..."


# ------------------------------------------------------------
# Markdown report
# ------------------------------------------------------------

def mutation_table_to_markdown(mutations_df):
    if mutations_df.empty:
        return "No mutations detected.\n"

    lines = []
    lines.append(
        "| ID | Type | Reference position | Reference sequence | Sample sequence | Length | Note |"
    )
    lines.append("|---:|---|---|---|---|---:|---|")

    for _, row in mutations_df.iterrows():
        lines.append(
            "| "
            f"{row['mutation_id']} | "
            f"{row['event_type']} | "
            f"{row['reference_position']} | "
            f"{shorten_text(row['reference_sequence'], 25)} | "
            f"{shorten_text(row['sample_sequence'], 25)} | "
            f"{row['length']} | "
            f"{shorten_text(row['note'], 80)} |"
        )

    return "\n".join(lines) + "\n"


def features_table_to_markdown(features):
    if not features:
        return "No selected GenBank features found.\n"

    lines = []
    lines.append("| Type | Label | Start | End | Strand |")
    lines.append("|---|---|---:|---:|---|")

    for feature in features:
        lines.append(
            "| "
            f"{feature['type']} | "
            f"{feature['label']} | "
            f"{feature['start']} | "
            f"{feature['end']} | "
            f"{feature['strand']} |"
        )

    return "\n".join(lines) + "\n"


def create_markdown_report(
    record,
    mutations_df,
    topology,
    genbank_path,
    mutation_csv_path,
    alignment_path,
    map_image_path,
    output_md,
    project_name=None,
):
    plasmid_name = get_plasmid_name(record)
    plasmid_length = len(record)
    mutation_summary = summarize_mutations(mutations_df)
    features = get_top_features(record)

    report_date = datetime.now().strftime("%Y-%m-%d %H:%M")

    if not project_name:
        project_name = f"{plasmid_name} mutation QC report"

    alignment_status = (
        alignment_path if alignment_path and os.path.exists(alignment_path) else "Not included"
    )

    map_status = (
        map_image_path if map_image_path and os.path.exists(map_image_path) else "Not included"
    )

    md = f"""# {project_name}

Generated: {report_date}

## Overview

This report compares a reference plasmid sequence with a mutated or sample sequence and summarizes detected sequence differences.

This report is intended for research-use automation and visualization. Final biological interpretation should be reviewed by the researcher.

## Input files

| Item | Path |
|---|---|
| Reference GenBank | `{genbank_path}` |
| Mutation report CSV | `{mutation_csv_path}` |
| Alignment preview | `{alignment_status}` |
| Plasmid map image | `{map_status}` |

## Plasmid information

| Field | Value |
|---|---|
| Plasmid name | {plasmid_name} |
| Length | {plasmid_length:,} bp |
| Topology | {topology} |
| GenBank features | {len(record.features)} |

## Mutation summary

| Type | Count |
|---|---:|
| Total mutations | {mutation_summary["total"]} |
| Mismatches | {mutation_summary["mismatch"]} |
| Insertions | {mutation_summary["insertion"]} |
| Deletions | {mutation_summary["deletion"]} |
| Other | {mutation_summary["other"]} |

## Detected mutations

{mutation_table_to_markdown(mutations_df)}

## Selected GenBank features

{features_table_to_markdown(features)}

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
"""

    output_dir = os.path.dirname(output_md)

    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    with open(output_md, "w", encoding="utf-8") as file:
        file.write(md)

    return output_md


# ------------------------------------------------------------
# PDF report
# ------------------------------------------------------------

def create_pdf_report(
    record,
    mutations_df,
    topology,
    genbank_path,
    mutation_csv_path,
    alignment_path,
    map_image_path,
    output_pdf,
    project_name=None,
):
    try:
        from reportlab.lib import colors
        from reportlab.lib.pagesizes import letter
        from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
        from reportlab.lib.units import inch
        from reportlab.platypus import (
            SimpleDocTemplate,
            Paragraph,
            Spacer,
            Table,
            TableStyle,
            Image,
            PageBreak,
            Preformatted,
        )
    except ImportError:
        raise ImportError(
            "reportlab is not installed. Run: pip install reportlab"
        )

    output_dir = os.path.dirname(output_pdf)

    if output_dir:
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    plasmid_name = get_plasmid_name(record)
    plasmid_length = len(record)
    mutation_summary = summarize_mutations(mutations_df)
    features = get_top_features(record)
    report_date = datetime.now().strftime("%Y-%m-%d %H:%M")

    if not project_name:
        project_name = f"{plasmid_name} mutation QC report"

    doc = SimpleDocTemplate(
        output_pdf,
        pagesize=letter,
        rightMargin=0.6 * inch,
        leftMargin=0.6 * inch,
        topMargin=0.6 * inch,
        bottomMargin=0.6 * inch,
    )

    styles = getSampleStyleSheet()

    styles.add(
        ParagraphStyle(
            name="SmallText",
            parent=styles["Normal"],
            fontSize=8,
            leading=10,
        )
    )

    story = []

    story.append(Paragraph(project_name, styles["Title"]))
    story.append(Spacer(1, 0.15 * inch))
    story.append(Paragraph(f"Generated: {report_date}", styles["Normal"]))
    story.append(Spacer(1, 0.2 * inch))

    story.append(Paragraph("Overview", styles["Heading2"]))
    story.append(
        Paragraph(
            "This report summarizes detected sequence differences between a "
            "reference plasmid and a mutated/sample sequence. It is intended "
            "for research-use automation and visualization. Final biological "
            "interpretation should be reviewed by the researcher.",
            styles["Normal"],
        )
    )
    story.append(Spacer(1, 0.2 * inch))

    story.append(Paragraph("Plasmid information", styles["Heading2"]))

    plasmid_info_table = Table(
        [
            ["Field", "Value"],
            ["Plasmid name", plasmid_name],
            ["Length", f"{plasmid_length:,} bp"],
            ["Topology", topology],
            ["GenBank features", str(len(record.features))],
            ["Reference GenBank", genbank_path],
            ["Mutation CSV", mutation_csv_path],
        ],
        colWidths=[1.8 * inch, 5.0 * inch],
    )

    plasmid_info_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ("FONTSIZE", (0, 0), (-1, -1), 8),
            ]
        )
    )

    story.append(plasmid_info_table)
    story.append(Spacer(1, 0.25 * inch))

    story.append(Paragraph("Mutation summary", styles["Heading2"]))

    summary_table = Table(
        [
            ["Type", "Count"],
            ["Total mutations", mutation_summary["total"]],
            ["Mismatches", mutation_summary["mismatch"]],
            ["Insertions", mutation_summary["insertion"]],
            ["Deletions", mutation_summary["deletion"]],
            ["Other", mutation_summary["other"]],
        ],
        colWidths=[2.5 * inch, 1.2 * inch],
    )

    summary_table.setStyle(
        TableStyle(
            [
                ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                ("GRID", (0, 0), (-1, -1), 0.4, colors.grey),
                ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                ("FONTSIZE", (0, 0), (-1, -1), 9),
            ]
        )
    )

    story.append(summary_table)
    story.append(Spacer(1, 0.25 * inch))

    if map_image_path and os.path.exists(map_image_path):
        story.append(Paragraph("Plasmid map", styles["Heading2"]))

        try:
            img = Image(map_image_path)
            img._restrictSize(6.8 * inch, 4.8 * inch)
            story.append(img)
            story.append(Spacer(1, 0.25 * inch))
        except Exception:
            story.append(
                Paragraph(
                    f"Map image could not be embedded: {map_image_path}",
                    styles["Normal"],
                )
            )

    story.append(Paragraph("Detected mutations", styles["Heading2"]))

    if mutations_df.empty:
        story.append(Paragraph("No mutations detected.", styles["Normal"]))
    else:
        mutation_rows = [
            [
                "ID",
                "Type",
                "Position",
                "Reference",
                "Sample",
                "Length",
            ]
        ]

        for _, row in mutations_df.iterrows():
            mutation_rows.append(
                [
                    str(row["mutation_id"]),
                    str(row["event_type"]),
                    str(row["reference_position"]),
                    shorten_text(row["reference_sequence"], 18),
                    shorten_text(row["sample_sequence"], 18),
                    str(row["length"]),
                ]
            )

        mutation_table = Table(
            mutation_rows,
            repeatRows=1,
            colWidths=[
                0.45 * inch,
                0.85 * inch,
                1.1 * inch,
                1.4 * inch,
                1.4 * inch,
                0.6 * inch,
            ],
        )

        mutation_table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                    ("GRID", (0, 0), (-1, -1), 0.35, colors.grey),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, -1), 7),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ]
            )
        )

        story.append(mutation_table)

    story.append(Spacer(1, 0.25 * inch))

    story.append(Paragraph("Selected GenBank features", styles["Heading2"]))

    if not features:
        story.append(
            Paragraph("No selected GenBank features found.", styles["Normal"])
        )
    else:
        feature_rows = [["Type", "Label", "Start", "End", "Strand"]]

        for feature in features:
            feature_rows.append(
                [
                    feature["type"],
                    shorten_text(feature["label"], 30),
                    str(feature["start"]),
                    str(feature["end"]),
                    str(feature["strand"]),
                ]
            )

        feature_table = Table(
            feature_rows,
            repeatRows=1,
            colWidths=[
                1.1 * inch,
                2.5 * inch,
                0.8 * inch,
                0.8 * inch,
                0.6 * inch,
            ],
        )

        feature_table.setStyle(
            TableStyle(
                [
                    ("BACKGROUND", (0, 0), (-1, 0), colors.lightgrey),
                    ("GRID", (0, 0), (-1, -1), 0.35, colors.grey),
                    ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
                    ("FONTSIZE", (0, 0), (-1, -1), 7),
                    ("VALIGN", (0, 0), (-1, -1), "TOP"),
                ]
            )
        )

        story.append(feature_table)

    alignment_preview = read_alignment_preview(alignment_path)

    if alignment_preview:
        story.append(PageBreak())
        story.append(Paragraph("Alignment preview", styles["Heading2"]))
        story.append(
            Preformatted(
                alignment_preview,
                styles["Code"],
            )
        )

    story.append(Spacer(1, 0.25 * inch))
    story.append(Paragraph("Limitations", styles["Heading2"]))
    story.append(
        Paragraph(
            "- Mutation positions are reported using the original/reference "
            "plasmid coordinate system.<br/>"
            "- Insertions and deletions in repetitive regions may have ambiguous "
            "exact positions depending on alignment scoring.<br/>"
            "- This script does not provide clinical or diagnostic interpretation.<br/>"
            "- Large rearrangements may require manual review or a more advanced "
            "alignment workflow.",
            styles["SmallText"],
        )
    )

    doc.build(story)

    return output_pdf


# ------------------------------------------------------------
# Main pipeline
# ------------------------------------------------------------

def generate_summary_report(
    genbank_path,
    mutation_csv_path,
    alignment_path,
    map_image_path,
    output_md,
    output_pdf,
    topology,
    project_name=None,
):
    print("Loading input files...")

    record = read_genbank_record(genbank_path)
    mutations_df = read_mutation_report(mutation_csv_path)

    print(f"Reference plasmid: {get_plasmid_name(record)}")
    print(f"Reference length: {len(record):,} bp")
    print(f"Mutations found: {len(mutations_df)}")

    print("\nGenerating Markdown report...")
    create_markdown_report(
        record=record,
        mutations_df=mutations_df,
        topology=topology,
        genbank_path=genbank_path,
        mutation_csv_path=mutation_csv_path,
        alignment_path=alignment_path,
        map_image_path=map_image_path,
        output_md=output_md,
        project_name=project_name,
    )
    print(f"✓ Markdown report saved: {output_md}")

    print("\nGenerating PDF report...")
    create_pdf_report(
        record=record,
        mutations_df=mutations_df,
        topology=topology,
        genbank_path=genbank_path,
        mutation_csv_path=mutation_csv_path,
        alignment_path=alignment_path,
        map_image_path=map_image_path,
        output_pdf=output_pdf,
        project_name=project_name,
    )
    print(f"✓ PDF report saved: {output_pdf}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate Markdown and PDF plasmid mutation QC summary report."
    )

    parser.add_argument(
        "--genbank",
        default="data/reference/reference_puc19.gb",
        help="Reference GenBank file path",
    )

    parser.add_argument(
        "--mutations",
        default="outputs/tables/mutation_report.csv",
        help="Mutation report CSV path",
    )

    parser.add_argument(
        "--alignment",
        default="outputs/reports/alignment.txt",
        help="Optional alignment.txt path",
    )

    parser.add_argument(
        "--map-image",
        default="outputs/figures/puc19_mutation_map.png",
        help="Optional plasmid map PNG path",
    )

    parser.add_argument(
        "--output-md",
        default="outputs/reports/summary_report.md",
        help="Output Markdown report path",
    )

    parser.add_argument(
        "--output-pdf",
        default="outputs/reports/summary_report.pdf",
        help="Output PDF report path",
    )

    parser.add_argument(
        "--topology",
        default="circular",
        choices=["circular", "linear"],
        help="Plasmid topology",
    )

    parser.add_argument(
        "--project-name",
        default=None,
        help="Optional project/report title",
    )

    args = parser.parse_args()

    generate_summary_report(
        genbank_path=args.genbank,
        mutation_csv_path=args.mutations,
        alignment_path=args.alignment,
        map_image_path=args.map_image,
        output_md=args.output_md,
        output_pdf=args.output_pdf,
        topology=args.topology,
        project_name=args.project_name,
    )


if __name__ == "__main__":
    main()
    