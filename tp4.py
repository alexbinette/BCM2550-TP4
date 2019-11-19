"""Travail Pratique 4 - Alexandre Binette"""


import sys
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Importe le fichier au argv[1] et attrape les erreurs de filepath
try:
    file = sys.argv[1]
except (IOError, IndexError):
    error_message = """Error : Genbank file missing at argv[1],
    please specify filepath"""
    print(error_message, file=sys.stderr)
    exit

# Crée le CDS avec les exons de la séquence Genbank
cds = ""

with open(file, "r") as f:
    for record in SeqIO.parse(f, "genbank"):
        for feature in record.features:
            if feature.type == "exon":
                cds += record.seq[feature.location.start: feature.location.end]


# Utilise un regex pour trouver les codons START de la séquence protéique
longest_prot = ""

for match in re.finditer('ATG', str(cds)):
    cds_trans = cds[match.start():].translate(to_stop=True)
    if len(longest_prot) < len(cds_trans):
        longest_prot = (cds_trans)

# Écrit le SeqRecord avec les paramètres nécessaires
prot_id = "NP_001107.2"
prot_desc = "adcy9_protein.fasta"

seq = SeqRecord(longest_prot, id = prot_id, description = prot_desc)

# Écrit le fichier FASTA
with open("./adcy9_protein.fasta", "w") as f:
    SeqIO.write(seq, f, format="fasta")
