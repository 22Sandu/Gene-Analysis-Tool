#!/usr/bin/env python
# standard library
from __future__ import print_function
import sys
import json
import os
import uuid

from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML

# program = sys.argv[1]
# database = sys.argv[2]
sequence = sys.argv[3]
# outputType = sys.argv[4]

asnFileName = str(uuid.uuid4()) + '.asn'
xmlFileName = str(uuid.uuid4()) + '.xml'
textFileName = str(uuid.uuid4()) + '.text'
outFileName = str(uuid.uuid4()) + '.fa'
inFileName = str(uuid.uuid4()) + '.fa'

inFile = open(inFileName, 'w+')
inFile.write(sequence)
inFile.close()

queryFileName = inFileName

dbPath = '/Users/anuradhawick/Desktop/blastpdb/sorted_env_nr'

os.system(
    'blastp -query ' + queryFileName + ' -num_threads 8 -max_hsps 10  -db ' + dbPath + ' -outfmt 11 > ' + asnFileName)
os.system('blast_formatter -archive ' + asnFileName + ' -outfmt 5 -out ' + xmlFileName)
os.system('blast_formatter -archive ' + asnFileName + ' -outfmt 0 -out ' + textFileName)

result_handle = open(xmlFileName)
blast_record = NCBIXML.read(result_handle)

outFile = open(outFileName, 'w+')

for alignment in blast_record.alignments:
    for hsp in alignment.hsps:
        outFile.write('>' + alignment.title + '\n')
        outFile.write(hsp.sbjct.replace("-", "") + '\n')
outFile.close()

outFile = open(outFileName)
out = outFile.read()
outFile.close()

outTextFile = open(textFileName)
outText = outTextFile.read()
outTextFile.close()

os.system('rm ' + xmlFileName)
os.system('rm ' + outFileName)
os.system('rm ' + asnFileName)
os.system('rm ' + textFileName)
os.system('rm ' + inFileName)

# blast_results = result_handle.read()

# print((blast_results))
print(json.dumps({'text': outText, 'output': out}))
