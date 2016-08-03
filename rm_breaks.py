#!/usr/bin/env python3
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-f", "--fasta", dest="fasta",
                    help="Fasta file with temporary names")

args = parser.parse_args()

with open(args.fasta, 'r') as fasta_records:
    for line in fasta_records:
        if '>' not in line:
            print(line.rstrip())
        else:
            print(line)
