#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fasta", dest="fasta",
                    help="Fasta file with temporary names")
parser.add_argument("-n", "--names", dest="names",
                    help="Name replacements for fasta file")

args = parser.parse_args()

replacement_names = open(args.names, 'r')
with open(args.fasta, 'r') as fasta_records:
    for line in fasta_records:
        if ">" in line:
            flag = False
            for name in replacement_names:
                gene_name = ' '.join(name.split('\t')[1:]).rstrip()
                gene_id = name.split('\t')[0].rstrip()
                if gene_id in line:
                    print(">" + name.rstrip())
                    flag = True
                    break
            if flag == False:
                print("ERROR: GeneID not found in names file" + line.rstrip())
        else:
            print(line.rstrip())
replacement_names.close()
