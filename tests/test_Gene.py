import pytest
import pickle
import numpy as np
import pandas as pd

from sequences import Gene


with open("data/BW25113_genome.pkl", "rb") as input:
    GENOME = pickle.load(input)

CODON_TABLE = pd.read_csv("data/Edinburgh-Genome-Foundry_e_coli_316407.csv")

# Fixtures for each test:
@pytest.fixture
def gene_fwd():
    return Gene(name="gene1", sequence="ATGAGATAA", strand=1)

@pytest.fixture
def gene_rev():
    return Gene(name="gene2", sequence="TTATGGCAT", strand=-1)

@pytest.fixture
def incomplete_fwd():
    return Gene(name="inc_gene1", sequence="ATGAGATA", strand=1)

@pytest.fixture
def no_gene():
    return Gene(name="no_gene", sequence="", strand=1)

@pytest.fixture
def gene_fwd_optimal():
    return Gene(name="gene3", sequence="ATGTGCTAA", strand=1)


# Test reverse complement sequence:
def test_revcompl_fwd(gene_fwd):
    assert gene_fwd.reverse_complement() == "TTATCTCAT"

def test_revcompl_rev(gene_rev):
    assert gene_rev.reverse_complement() == "ATGCCATAA"

def test_revcompl_inc(incomplete_fwd):
    assert incomplete_fwd.reverse_complement() == "TATCTCAT"

def test_revcompl_no(no_gene):
    assert no_gene.reverse_complement() == ""


# Test transcription:
def test_transcribe_fwd(gene_fwd):
    assert gene_fwd.transcribe() == "AUGAGAUAA"

def test_transcribe_rev(gene_rev):
    assert gene_rev.transcribe() == "AUGCCAUAA"

def test_transcribe_inc(incomplete_fwd):
    assert incomplete_fwd.transcribe() == "AUGAGAUA"

def test_transcribe_no(no_gene):
    assert no_gene.transcribe() == ""


# Test transcribed sequence
def test_transcribedseq_fwd(gene_fwd):
    assert gene_fwd.transcribed_sequence == "AUGAGAUAA"

def test_transcribedseq_rev(gene_rev):
    assert gene_rev.transcribed_sequence == "AUGCCAUAA"

def test_transcribedseq_inc(incomplete_fwd):
    assert incomplete_fwd.transcribed_sequence == "AUGAGAUA"

def test_transcribeseq_no(no_gene):
    assert no_gene.transcribed_sequence == ""


## Test GC content calculation:
def test_gc_fwd(gene_fwd):
    assert gene_fwd.gc_content() == 2/9

def test_gc_rev(gene_rev):
    assert gene_rev.gc_content() == 3/9

def test_gc_inc(incomplete_fwd):
    assert incomplete_fwd.gc_content() == 2/8

def test_transcribeseq_no(no_gene):
    assert np.isnan(no_gene.gc_content())    


## Test no. codon check:
def test_codoncheck_fwd(gene_fwd):
    assert gene_fwd._check_sequence_length() is True

def test_codoncheck_rev(gene_rev):
    assert gene_rev._check_sequence_length() is True

def test_codoncheck_inc(incomplete_fwd):
    assert incomplete_fwd._check_sequence_length() is False

def test_codoncheck_no(no_gene):
    assert no_gene._check_sequence_length() is True


## Test no. codon determination:
def test_nocodons_fwd(gene_fwd):
    assert gene_fwd.no_codons == 3

def test_nocodons_rev(gene_rev):
    assert gene_rev.no_codons == 3

def test_nocodons_inc(incomplete_fwd):
    assert np.isnan(incomplete_fwd.no_codons)

def test_nocodons_no(no_gene):
    assert no_gene.no_codons == 0


## Test codon determination:
def test_codonseq_fwd(gene_fwd):
    assert gene_fwd.codon_sequence == ["AUG", "AGA", "UAA"]

def test_codonseq_rev(gene_rev):
    assert gene_rev.codon_sequence == ["AUG", "CCA", "UAA"]

def test_codonseq_inc(incomplete_fwd):
    assert incomplete_fwd.codon_sequence == []

def test_codonseq_no(no_gene):
    assert no_gene.codon_sequence == []


## Test codon determination:
def test_cai_fwd(gene_fwd):
    assert np.isclose(gene_fwd.codon_adaptation_index(CODON_TABLE), 0.464159)

def test_cai_fwd_opt(gene_fwd_optimal):
    assert np.isclose(gene_fwd_optimal.codon_adaptation_index(CODON_TABLE), 1.00)

def test_cai_rev(gene_rev):
    assert np.isclose(gene_rev.codon_adaptation_index(CODON_TABLE), 0.710383)

def test_cai_inc(incomplete_fwd):
    assert np.isnan(incomplete_fwd.codon_adaptation_index(CODON_TABLE))

def test_cai_no(no_gene):
    assert np.isnan(no_gene.codon_adaptation_index(CODON_TABLE))