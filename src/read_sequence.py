#goal: read FASTA or GENBANK file and return the sequence

from Bio import SeqIO

def read_fasta_sequence(file_path):
    record = SeqIO.read(file_path, "fasta")
    return str(record.seq)

def read_genbank_record(file_path):
    record = SeqIO.read(file_path, "genbank")
    return record