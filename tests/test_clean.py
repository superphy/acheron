import pytest
import os

from acheron.workflows import clean
test_fasta = 'data/acheron_test_samples/wgs/raw/subSRR2407535.fasta'
out = 'data/acheron_test_samples/wgs/clean/'

def test_get_files_to_analyze():
    assert False
    files = clean.get_files_to_analyze(test_fasta)
    assert len(files) == 1
    assert '/'.join(files[0].split('/')[-5:]) == test_fasta

def test_find_recurring_char():
    record = 'AAAT'
    assert clean.find_recurring_char(0.70, 'AAAT', 0, 4) == 'A'
    assert clean.find_recurring_char(0.80, 'AAAT', 0, 4) == 'X'

def test_format_files():
    files = clean.get_files_to_analyze(test_fasta)
    os.makedirs(out, exist_ok=True)
    clean.format_files(files, out)
    from Bio import SeqIO

    num_contigs = 0
    with open(files[0].replace("raw","clean"), 'r') as fh:
        for _ in SeqIO.parse(fh, "fasta"):
            num_contigs+=1
    assert num_contigs == 3
