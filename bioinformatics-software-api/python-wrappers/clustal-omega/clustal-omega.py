#!/usr/bin/env python
# standard library
from Bio.Align.Applications import ClustalOmegaCommandline

#Input file
in_file = "Acanthaster_planci_Gnomon.fasta"

#Output file
out_file = "aligned.fasta"

#Get Clustal omega command
clustalomega_cline = ClustalOmegaCommandline(infile=in_file, outfile=out_file, verbose=True, auto=True)
print(clustalomega_cline)
