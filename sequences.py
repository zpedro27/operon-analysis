import collections
import numpy as np
import pandas as pd
from typing import List

BASE_PAIRING = {"A": "T",
                "G": "C", 
                "C": "G",
                "T": "A"}

def observed_freq(codon: str, codon_table: pd.DataFrame):
    try:
        codon_frequency, = codon_table.loc[codon_table["codon"] == codon, "relative_frequency"]
    except ValueError:
        codon_frequency = np.nan
    return codon_frequency

def max_freq(codon: str, codon_table: pd.DataFrame):
    try:
        coded_aminoacid, = codon_table.loc[codon_table["codon"] == codon, "amino_acid"].unique()
        max_codon_frequency = codon_table.loc[codon_table["amino_acid"] == coded_aminoacid, "relative_frequency"].max()
    except ValueError:
        max_codon_frequency = np.nan
    return max_codon_frequency

class GenomeSequence:
    def __init__(self, name: str, sequence: str, strand: int) -> None:
        self._name = name
        self._sequence = sequence
        self._strand = strand
        self._transcribed_sequence = self.transcribe()
    
    def __str__(self):
        return f"{self.__class__.__name__}: {self.name}, {self.length} bp, strand {self.strand}"

    @property
    def name(self):
        return self._name

    @property
    def sequence(self):
        return self._sequence
    
    @property
    def transcribed_sequence(self):
        return self._transcribed_sequence
    
    @property
    def strand(self):
        return self._strand
    
    @property
    def length(self):
        return len(self.sequence)

    def gc_content(self):
        base_counter = collections.Counter(self.sequence)
        total_gc = base_counter["G"] + base_counter["C"]
        if self.length > 0:
            return total_gc / self.length
        return np.nan

    def reverse_complement(self):
        paired = [BASE_PAIRING[base] for base in self.sequence]
        rev_paired  = paired[::-1]
        return "".join(rev_paired)

    def transcribe(self):
        if self.strand == -1:
            return self.reverse_complement().replace("T", "U")
        return self.sequence.replace("T", "U")


class Gene(GenomeSequence):
    def __init__(self, name: str, sequence: str, strand: int) -> None:
        super().__init__(name=name, sequence=sequence, strand=strand)
        self._divisible_by_three = self._check_sequence_length()

    def _check_sequence_length(self):
        if self.length % 3 != 0:
            print(f"Gene sequence ({self.name}) is not divisible by 3!")
            return False
        return True

    @property
    def no_codons(self):
        if self._divisible_by_three is False:
            return np.nan
        return int(self.length/3)
        

    @property
    def codon_sequence(self):
        if self._divisible_by_three is False:
            return []
        return [self.transcribed_sequence[3*i : 3*(i+1)] for i in range(self.no_codons)]
        

    def codon_adaptation_index(self, codon_table: pd.DataFrame):
        relative_adaptiveness = [ observed_freq(codon, codon_table) / max_freq(codon, codon_table) for codon in self.codon_sequence]
        if (self._divisible_by_three is False) or (self.no_codons == 0):
            return np.nan
        return np.prod(relative_adaptiveness)**(1/self.no_codons)
    

class Operon(GenomeSequence):
    def __init__(self, name: str, sequence: str, strand: int, genes: List[Gene]) -> None:
        super().__init__(name=name, sequence=sequence, strand=strand)
        self._genes = genes
    
    @property
    def genes(self):
        return self._genes

    def spacer_length(self):
        total_length = self.length
        total_orf_length = sum([gene.length for gene in self.genes])
        return total_length - total_orf_length
    
    def codon_adaptation_index(self, codon_table: dict):
        return {gene.name: gene.codon_adaptation_index(codon_table) for gene in self.genes}



