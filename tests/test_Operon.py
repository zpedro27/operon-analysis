import pytest
import pickle
import numpy as np

from sequences import Gene, Operon

PATH_TO_GENOME = r".................\\Data\\General\\"
with open(PATH_TO_GENOME+"BW25113_genome.pkl", "rb") as input:
    GENOME = pickle.load(input)

# Fixtures for each test:
@pytest.fixture
def operon_fwd():
    seq1 = "ATGAGATAA"
    seq2 = "ATGGCGTAA"
    seq_operon = f"{seq1}GGAC{seq2}"
    gene1 = Gene(name="gene1", sequence=seq1, strand=1)
    gene2 = Gene(name="gene2", sequence=seq2, strand=1)
    return Operon(name="operon1", sequence=seq_operon, genes=[gene1, gene2], strand=1)

@pytest.fixture
def operon_rev():
    seq1 = "TTATCTCAT"
    seq2 = "TTAGGGCAT"
    seq_operon = f"{seq1}CAGG{seq2}"
    gene1 = Gene(name="gene1", sequence=seq1, strand=-1)
    gene2 = Gene(name="gene2", sequence=seq2, strand=-1)
    return Operon(name="operon2", sequence=seq_operon, genes=[gene1, gene2], strand=-1)

@pytest.fixture
def operon_nospacer_fwd():
    seq1 = "ATGAGATAA"
    seq2 = "ATGGCGTAA"
    seq_operon = f"{seq1}{seq2}"
    gene1 = Gene(name="gene1", sequence=seq1, strand=1)
    gene2 = Gene(name="gene2", sequence=seq2, strand=1)
    return Operon(name="operon3", sequence=seq_operon, genes=[gene1, gene2], strand=1)


# Test reverse complement sequence:
def test_revcompl_fwd(operon_fwd):
    assert operon_fwd.reverse_complement() == "TTACGCCATGTCCTTATCTCAT"

def test_revcompl_nospacer(operon_nospacer_fwd):
    assert operon_nospacer_fwd.reverse_complement() == "TTACGCCATTTATCTCAT"

def test_revcompl_rev(operon_rev):
    assert operon_rev.reverse_complement() == "ATGCCCTAACCTGATGAGATAA"


# Test spacer length determination
def test_spacer_fwd(operon_fwd):
    assert operon_fwd.spacer_length() == 4

def test_spacer_nospacer(operon_nospacer_fwd):
    assert operon_nospacer_fwd.spacer_length() == 0

def test_spacer_rev(operon_rev):
    assert operon_rev.spacer_length() == 4
