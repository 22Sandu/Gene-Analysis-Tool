from Bio.Align.Applications import TCoffeeCommandline
tcoffee_cline = TCoffeeCommandline(infile="aligned.fasta",
                                   output="clustalw",
                                   outfile="aligned.aln")
print(tcoffee_cline)

#in_file = "/Users/vijinimallawaarachchi/Documents/Python/Acanthaster_planci_Gnomon.fsa"

#t_coffee -output clustalw -infile aligned.fasta -outfile aligned.aln
