"""
goal: take original sequence and create mutated sequence
output: data/mutated/mutated_puc19.fasta
"""

#!/usr/bin/env python3
"""
Simulate mutations on DNA sequences.
Handles substitutions, deletions, and insertions with proper
conversion between biological (1-based) and Python (0-based) positions.
"""

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml
import os


def read_fasta_sequence(file_path):
    """
    Read one FASTA file and return sequence string.
    """
    record = SeqIO.read(file_path, "fasta")
    return str(record.seq).upper(), record.id, record.description


def write_fasta_sequence(sequence, output_path, record_id="mutated_sequence", description="Simulated mutated plasmid"):
    """
    Save sequence as FASTA file.
    """
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    record = SeqRecord(
        Seq(sequence),
        id=record_id,
        description=description
    )

    SeqIO.write(record, output_path, "fasta")


def get_mutation_sort_position(mutation):
    """
    Get position for sorting mutations.
    This is important because insertions/deletions change sequence length.
    """
    mutation_type = mutation["type"]

    if mutation_type == "change":
        return mutation["position"]

    if mutation_type == "delete":
        return mutation["start"]

    if mutation_type == "insert":
        return mutation["after"]

    raise ValueError(f"Unknown mutation type: {mutation_type}")


def simulate_mutations(original_sequence, mutations):
    """
    Apply mutations to original sequence.

    Important:
    - Positions in config.yml are biological positions.
    - Biological positions start from 1.
    - Python string indexes start from 0.
    - Mutation positions are based on the ORIGINAL plasmid sequence.
    """

    mutated_sequence = original_sequence
    mutation_log = []

    # Sort mutations by original plasmid position
    mutations = sorted(mutations, key=get_mutation_sort_position)

    # Offset tracks sequence length changes after insertions/deletions
    offset = 0

    for mutation in mutations:
        mutation_type = mutation["type"]

        if mutation_type == "change":
            position = mutation["position"]
            new_base = mutation["new_base"].upper()

            index = position - 1 + offset

            old_base = mutated_sequence[index]

            mutated_sequence = (
                mutated_sequence[:index]
                + new_base
                + mutated_sequence[index + 1:]
            )

            mutation_log.append({
                "type": "change",
                "position": position,
                "old_base": old_base,
                "new_base": new_base,
                "note": f"{old_base} changed to {new_base} at position {position}"
            })

        elif mutation_type == "delete":
            start = mutation["start"]
            end = mutation["end"]

            start_index = start - 1 + offset
            end_index = end + offset

            deleted_sequence = mutated_sequence[start_index:end_index]

            mutated_sequence = (
                mutated_sequence[:start_index]
                + mutated_sequence[end_index:]
            )

            deleted_length = end - start + 1
            offset -= deleted_length

            mutation_log.append({
                "type": "delete",
                "start": start,
                "end": end,
                "deleted_sequence": deleted_sequence,
                "length": deleted_length,
                "note": f"{deleted_length} bp deleted from position {start} to {end}"
            })

        elif mutation_type == "insert":
            after = mutation["after"]
            insert_sequence = mutation["sequence"].upper()

            insert_index = after + offset

            mutated_sequence = (
                mutated_sequence[:insert_index]
                + insert_sequence
                + mutated_sequence[insert_index:]
            )

            offset += len(insert_sequence)

            mutation_log.append({
                "type": "insert",
                "after": after,
                "inserted_sequence": insert_sequence,
                "length": len(insert_sequence),
                "note": f"{len(insert_sequence)} bp inserted after position {after}"
            })

        else:
            raise ValueError(f"Unknown mutation type: {mutation_type}")

    return mutated_sequence, mutation_log


def load_config(config_path):
    """
    Load config.yml file.
    """
    with open(config_path, "r") as file:
        return yaml.safe_load(file)


def run_simulation(config_path="config.yml"):
    """
    Main function:
    1. Read config.yml
    2. Read original FASTA
    3. Apply mutations
    4. Save mutated FASTA
    5. Print mutation log
    """

    config = load_config(config_path)

    reference_fasta = config["reference_fasta"]
    output_fasta = config["mutated_fasta"]
    mutations = config["mutations"]

    original_sequence, record_id, description = read_fasta_sequence(reference_fasta)

    print(f"Original plasmid: {record_id}")
    print(f"Original length: {len(original_sequence)} bp")

    mutated_sequence, mutation_log = simulate_mutations(original_sequence, mutations)

    write_fasta_sequence(
        mutated_sequence,
        output_fasta,
        record_id=f"{record_id}_mutated",
        description="Simulated mutated plasmid sequence"
    )

    print(f"Mutated length: {len(mutated_sequence)} bp")
    print(f"Mutated FASTA saved to: {output_fasta}")
    print("\nMutation log:")

    for item in mutation_log:
        print("-", item["note"])

    return mutated_sequence, mutation_log


if __name__ == "__main__":
    run_simulation("config.yml")