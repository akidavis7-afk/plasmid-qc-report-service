from Bio import SeqIO
from Bio.Align import PairwiseAligner
import pandas as pd
import argparse
import os


def read_fasta_sequence(file_path):
    """
    Read a FASTA file and return:
    - sequence ID
    - sequence string
    """
    record = SeqIO.read(file_path, "fasta")
    sequence = str(record.seq).upper()

    return record.id, sequence


def build_gapped_alignment(reference_sequence, sample_sequence):
    """
    Align reference sequence and sample/mutated sequence.

    This creates two aligned strings with gaps:
    Example:

    Reference: ATGCC--TAA
    Sample:    ATG-CGGTA

    '-' means insertion or deletion.
    """

    aligner = PairwiseAligner()
    aligner.mode = "global"

    # Scoring settings
    # Match is good.
    # Mismatch is bad.
    # Gap is used for insertion/deletion.
    aligner.match_score = 2
    aligner.mismatch_score = -3
    aligner.open_gap_score = -5
    aligner.extend_gap_score = -1

    alignments = aligner.align(reference_sequence, sample_sequence)
    best_alignment = alignments[0]

    reference_blocks, sample_blocks = best_alignment.aligned

    gapped_reference = []
    gapped_sample = []

    reference_pos = 0
    sample_pos = 0

    for ref_block, sample_block in zip(reference_blocks, sample_blocks):
        ref_start, ref_end = ref_block
        sample_start, sample_end = sample_block

        # Reference has extra bases before this aligned block.
        # This means deletion in sample.
        if ref_start > reference_pos:
            deleted_part = reference_sequence[reference_pos:ref_start]
            gapped_reference.append(deleted_part)
            gapped_sample.append("-" * len(deleted_part))
            reference_pos = ref_start

        # Sample has extra bases before this aligned block.
        # This means insertion in sample.
        if sample_start > sample_pos:
            inserted_part = sample_sequence[sample_pos:sample_start]
            gapped_reference.append("-" * len(inserted_part))
            gapped_sample.append(inserted_part)
            sample_pos = sample_start

        # Aligned region
        ref_part = reference_sequence[ref_start:ref_end]
        sample_part = sample_sequence[sample_start:sample_end]

        gapped_reference.append(ref_part)
        gapped_sample.append(sample_part)

        reference_pos = ref_end
        sample_pos = sample_end

    # Remaining reference tail = deletion
    if reference_pos < len(reference_sequence):
        deleted_part = reference_sequence[reference_pos:]
        gapped_reference.append(deleted_part)
        gapped_sample.append("-" * len(deleted_part))

    # Remaining sample tail = insertion
    if sample_pos < len(sample_sequence):
        inserted_part = sample_sequence[sample_pos:]
        gapped_reference.append("-" * len(inserted_part))
        gapped_sample.append(inserted_part)

    return "".join(gapped_reference), "".join(gapped_sample)


def format_reference_position(event_type, start, end=None):
    """
    Format position based on original reference sequence coordinate.
    """
    if event_type == "insertion":
        if start == 0:
            return "before position 1"
        return f"after {start}"

    if start == end:
        return str(start)

    return f"{start}-{end}"


def detect_mutations(gapped_reference, gapped_sample):
    """
    Detect mutation events from gapped alignment.

    Reports positions based on the ORIGINAL reference sequence.
    """

    mutations = []
    reference_position = 0
    i = 0
    mutation_id = 1

    alignment_length = len(gapped_reference)

    while i < alignment_length:
        ref_base = gapped_reference[i]
        sample_base = gapped_sample[i]

        # Case 1: same base, no mutation
        if ref_base == sample_base:
            if ref_base != "-":
                reference_position += 1
            i += 1
            continue

        # Case 2: mismatch / base substitution
        if ref_base != "-" and sample_base != "-":
            start_position = reference_position + 1
            ref_seq = []
            sample_seq = []

            while (
                i < alignment_length
                and gapped_reference[i] != "-"
                and gapped_sample[i] != "-"
                and gapped_reference[i] != gapped_sample[i]
            ):
                ref_seq.append(gapped_reference[i])
                sample_seq.append(gapped_sample[i])
                reference_position += 1
                i += 1

            end_position = reference_position

            ref_seq = "".join(ref_seq)
            sample_seq = "".join(sample_seq)

            mutations.append({
                "mutation_id": mutation_id,
                "event_type": "mismatch",
                "reference_position": format_reference_position(
                    "mismatch",
                    start_position,
                    end_position
                ),
                "reference_start": start_position,
                "reference_end": end_position,
                "reference_sequence": ref_seq,
                "sample_sequence": sample_seq,
                "length": len(ref_seq),
                "note": f"{ref_seq} changed to {sample_seq} at reference position {format_reference_position('mismatch', start_position, end_position)}"
            })

            mutation_id += 1
            continue

        # Case 3: deletion
        # Reference has bases, sample has gaps.
        if ref_base != "-" and sample_base == "-":
            start_position = reference_position + 1
            deleted_seq = []

            while (
                i < alignment_length
                and gapped_reference[i] != "-"
                and gapped_sample[i] == "-"
            ):
                deleted_seq.append(gapped_reference[i])
                reference_position += 1
                i += 1

            end_position = reference_position
            deleted_seq = "".join(deleted_seq)

            mutations.append({
                "mutation_id": mutation_id,
                "event_type": "deletion",
                "reference_position": format_reference_position(
                    "deletion",
                    start_position,
                    end_position
                ),
                "reference_start": start_position,
                "reference_end": end_position,
                "reference_sequence": deleted_seq,
                "sample_sequence": "-",
                "length": len(deleted_seq),
                "note": f"{len(deleted_seq)} bp deleted from reference position {format_reference_position('deletion', start_position, end_position)}"
            })

            mutation_id += 1
            continue

        # Case 4: insertion
        # Reference has gap, sample has extra bases.
        if ref_base == "-" and sample_base != "-":
            inserted_after_position = reference_position
            inserted_seq = []

            while (
                i < alignment_length
                and gapped_reference[i] == "-"
                and gapped_sample[i] != "-"
            ):
                inserted_seq.append(gapped_sample[i])
                i += 1

            inserted_seq = "".join(inserted_seq)

            mutations.append({
                "mutation_id": mutation_id,
                "event_type": "insertion",
                "reference_position": format_reference_position(
                    "insertion",
                    inserted_after_position
                ),
                "reference_start": inserted_after_position,
                "reference_end": inserted_after_position,
                "reference_sequence": "-",
                "sample_sequence": inserted_seq,
                "length": len(inserted_seq),
                "note": f"{len(inserted_seq)} bp inserted {format_reference_position('insertion', inserted_after_position)}"
            })

            mutation_id += 1
            continue

    return mutations


def save_mutation_report(mutations, output_csv):
    """
    Save mutation report as CSV.
    """

    os.makedirs(os.path.dirname(output_csv), exist_ok=True)

    columns = [
        "mutation_id",
        "event_type",
        "reference_position",
        "reference_start",
        "reference_end",
        "reference_sequence",
        "sample_sequence",
        "length",
        "note"
    ]

    df = pd.DataFrame(mutations, columns=columns)
    df.to_csv(output_csv, index=False)

    return df


def save_alignment_file(gapped_reference, gapped_sample, output_file):
    """
    Save gapped alignment text file.
    This is useful for debugging and portfolio proof.
    """

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as file:
        file.write(">reference_aligned\n")
        file.write(gapped_reference + "\n")
        file.write(">sample_aligned\n")
        file.write(gapped_sample + "\n")


def compare_fasta_files(reference_fasta, sample_fasta, output_csv, alignment_output=None):
    """
    Main comparison function.
    """

    reference_id, reference_sequence = read_fasta_sequence(reference_fasta)
    sample_id, sample_sequence = read_fasta_sequence(sample_fasta)

    print(f"Reference file: {reference_fasta}")
    print(f"Reference ID: {reference_id}")
    print(f"Reference length: {len(reference_sequence)} bp")
    print()
    print(f"Sample file: {sample_fasta}")
    print(f"Sample ID: {sample_id}")
    print(f"Sample length: {len(sample_sequence)} bp")
    print()

    print("Running global alignment...")
    gapped_reference, gapped_sample = build_gapped_alignment(
        reference_sequence,
        sample_sequence
    )

    print("Detecting mutations...")
    mutations = detect_mutations(gapped_reference, gapped_sample)

    df = save_mutation_report(mutations, output_csv)

    if alignment_output:
        save_alignment_file(gapped_reference, gapped_sample, alignment_output)

    print(f"Detected mutations: {len(mutations)}")
    print(f"Mutation report saved to: {output_csv}")

    if alignment_output:
        print(f"Alignment file saved to: {alignment_output}")

    if len(df) > 0:
        print()
        print(df[[
            "mutation_id",
            "event_type",
            "reference_position",
            "reference_sequence",
            "sample_sequence",
            "length"
        ]].to_string(index=False))
    else:
        print("No mutations detected.")

    return df


def main():
    parser = argparse.ArgumentParser(
        description="Compare reference FASTA and mutated/sample FASTA to detect mismatches, insertions, and deletions."
    )

    parser.add_argument(
        "--reference",
        default="data/reference/reference_puc19.fasta",
        help="Path to reference FASTA file"
    )

    parser.add_argument(
        "--sample",
        default="data/mutated/mutated_puc19.fasta",
        help="Path to mutated/sample FASTA file"
    )

    parser.add_argument(
        "--output",
        default="outputs/tables/mutation_report.csv",
        help="Output CSV report path"
    )

    parser.add_argument(
        "--alignment-output",
        default="outputs/reports/alignment.txt",
        help="Optional output file for gapped alignment"
    )

    args = parser.parse_args()

    compare_fasta_files(
        reference_fasta=args.reference,
        sample_fasta=args.sample,
        output_csv=args.output,
        alignment_output=args.alignment_output
    )


if __name__ == "__main__":
    main()