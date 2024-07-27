import pytest
import pickle
from Bio.SeqFeature import FeatureLocation, CompoundLocation, ExactPosition

from locations import GenomeLocation

PATH_TO_GENOME = r"...............\\Data\\General\\"
with open(PATH_TO_GENOME+"BW25113_genome.pkl", "rb") as input:
    GENOME = pickle.load(input)

# Fixtures for each test:
@pytest.fixture
def single_fwd():
    location = [FeatureLocation(ExactPosition(350), ExactPosition(370), strand=1)]
    return GenomeLocation(name="seq1", location=location, genome=GENOME)

@pytest.fixture
def single_rev():
    location = [FeatureLocation(ExactPosition(9930), ExactPosition(9950), strand=-1)]
    return GenomeLocation(name="seq2", location=location, genome=GENOME)

@pytest.fixture
def list_fwd():
    location = [FeatureLocation(ExactPosition(350), ExactPosition(370), strand=1),
                FeatureLocation(ExactPosition(375), ExactPosition(380), strand=1)]
    return GenomeLocation(name="seq3", location=location, genome=GENOME)

@pytest.fixture
def compound_fwd():
    location = [CompoundLocation([FeatureLocation(ExactPosition(58470), ExactPosition(58480), strand=1),
                                 FeatureLocation(ExactPosition(58485), ExactPosition(58490), strand=1)],
                                "join")]
    return GenomeLocation(name="seq4", location=location, genome=GENOME)


# Tests:
def test_create_new_feature_singlefwd(single_fwd):
    to_compare = FeatureLocation(350, 370, strand=1)
    assert single_fwd._location == to_compare

def test_create_new_feature_singlerev(single_rev):
    to_compare = FeatureLocation(9930, 9950, strand=-1)
    assert single_rev._location == to_compare

def test_create_new_feature_listfwd(list_fwd):
    to_compare = FeatureLocation(350, 380, strand=1)
    assert list_fwd._location == to_compare

def test_create_new_feature_compfwd(compound_fwd):
    to_compare = FeatureLocation(58470, 58490, strand=1)
    assert compound_fwd._location == to_compare



def test_retrieve_seq_singlefwd(single_fwd):
    to_compare = "GTTCGGCGGTACATCAGTGG"
    assert single_fwd.retrieve_sequence() == to_compare

def test_retrieve_seq_singlerev(single_rev):
    to_compare = "GTGGGATTCACCAATCGGCA"   # not the reverse complement!
    assert single_rev.retrieve_sequence() == to_compare

def test_retrieve_seq_listfwd(list_fwd):
    to_compare = "GTTCGGCGGTACATCAGTGGCAAATGCAGA"
    assert list_fwd.retrieve_sequence() == to_compare

def test_retrieve_seq_compfwd(compound_fwd):
    to_compare = "ACCATGAAAGTATCAGTTCC"
    assert compound_fwd.retrieve_sequence() == to_compare
