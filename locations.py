from Bio.SeqFeature import FeatureLocation
from Bio.SeqRecord  import SeqRecord
from typing import List


class GenomeLocation:
    def __init__(self, name: str, location: List[FeatureLocation], genome: SeqRecord) -> None:
        self._name = name
        self._genome = genome
        self._location = self._create_new_feature(location)  #sum(location)

    @property
    def name(self):
        return self._name

    @property
    def start(self):
        return self._location.start
    
    @property
    def end(self):
        return self._location.end
    
    @property
    def strand(self):
        return self._location.strand
    
    def _create_new_feature(self, location: List[FeatureLocation]):
        whole_unit = sum(location)

        return FeatureLocation(whole_unit.start,
                               whole_unit.end,
                               strand=whole_unit.strand)

    def retrieve_sequence(self):
        return str(self._genome.seq[self.start : self.end])