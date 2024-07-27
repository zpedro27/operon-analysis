import pytest
import pandas as pd
import numpy as np

from sequences import observed_freq, max_freq


CODON_TABLE = pd.read_csv("data/Edinburgh-Genome-Foundry_e_coli_316407.csv")


# Test suite:
def test_observedfreq_aa():
    assert observed_freq("GAC", CODON_TABLE) == 0.37

def test_observedfreq_stop():
    assert observed_freq("UAA", CODON_TABLE) == 0.64

def test_observedfreq_start():
    assert observed_freq("AUG", CODON_TABLE) == 1.00

def test_observedfreq_none():
    assert np.isnan(observed_freq("", CODON_TABLE))


def test_maxfreq_aa():
    assert max_freq("GAC", CODON_TABLE) == 0.63

def test_maxfreq_stop():
    assert max_freq("UAA", CODON_TABLE) == 0.64

def test_maxfreq_start():
    assert max_freq("AUG", CODON_TABLE) == 1.00

def test_maxfreq_none():
    assert np.isnan(max_freq("", CODON_TABLE)) 