import pdb
import sys
from Bio import SeqIO

prev_geneID = ''
buffer = []

for record in SeqIO.parse(open(sys.argv[1]), "fasta"):
    geneID = str(record.description).split(' ')[3][0:]

    if prev_geneID == '':  # first pass
        prev_geneID = geneID
        buffer.append((str(record.description), str(record.seq)))
    else:
        if geneID != prev_geneID:
            buffer.sort(key=lambda x: len(x[1]))
            print('>' + buffer[-1][0])  # print gene info of top index
            print(buffer[-1][1])  # print sequence of top index
            prev_geneID = geneID
            buffer = [(str(record.description), str(record.seq))]
        else:
            buffer.append((str(record.description), str(record.seq)))

buffer.sort(key=lambda x: len(x[1]))
print('>' + buffer[-1][0])
print(buffer[-1][1])
