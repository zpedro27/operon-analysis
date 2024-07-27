from Bio.SeqRecord  import SeqRecord
from typing import List

def create_gene_location_table(genome: SeqRecord) -> dict:
    """
    Creates a dictionary of {"gene": location in genome}
    """
    hashtable = {}
    for feature in genome.features:
        if feature.type=="CDS":
            gene, = feature.qualifiers["gene"]
            hashtable[gene] = feature.location
    return hashtable 


def get_gene_names(unit_name: str, delimiter: str = ",") -> List[str]:
    """ 
    Given a string specifying a transcriptional unit, 
    determines the name of the first and last genes
    """
    return unit_name.split(delimiter)


def get_gene_instances(gene_names: List[str], table: dict) -> dict:
    """
    Retrieves a dictionary mapping gene names to their 
    corresponding SeqFeature objects.
    """
    return {gene: table[gene] for gene in gene_names}