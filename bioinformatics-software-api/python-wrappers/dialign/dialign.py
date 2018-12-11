from Bio.Align.Applications import DialignCommandline
dialign_cline = DialignCommandline(input="example.fasta", fn="aligned", fa=True)

print(dialign_cline)

#dialign2-2 -fa -fn aligned unaligned.fasta