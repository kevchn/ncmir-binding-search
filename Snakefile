#!/usr/bin/env python

"""
Author: K. Chen
Affiliation: NCI/NIH
Aim: A pipeline to find and analzye possible non-canonical miR targets using simple seed matching.
Date: Mon Nov  2 14:03:11 CET 2015
Run: snakemake
"""

import re
import json
import pandas as pd
import csv

###############################################################################

NAMES = ['cds','3utr','5utr']
MOTIF = "GCGGAAC"

###############################################################################

rule all:
    input:
        expand("output/{A}.found_genes.txt", A=NAMES),
        expand("summary_output/{A}.summary.txt", A=NAMES)

rule get_analysis_of_seqs:
    input:
        "genome/{A}.motif.fa"
    output:
        "output/{A}.found_genes.txt",
        "summary_output/{A}.summary.txt"
    run:
        try:    # Since there are multiple output files, deleting one before rerunning
                # Snakemake can cause the re-appending of output to the other file
            os.remove(output[0])
            os.remove(output[1])
        except OSError:
            pass
        with open(str(input), 'rt') as input_file:
            table_out = []
            prefix = ""
            num_genes_per_motif = {}
            for line in input_file:
                if '*****BEGIN' in line:
                    prefix = line.split(" ")[1].rstrip()
                    num_genes = 0  # reset
        # SECTION | GET INFO OF FOUND GENES ########################################
                if '>' in line:
                    num_genes += 1
                    gene_symbol = re.split(r'\t+', line)[0][1:]
                    try:
                        gene_desc = re.split(r'\t+', line)[1].rstrip()
                    except IndexError:
                        gene_desc = "No Description"
                    table_out.append([gene_symbol, gene_desc])

        # SECTION | DISPLAY FOUND GENES FOR EACH SITE ###############################
                if '*****END' in line:
                    num_genes_per_motif[prefix] = num_genes
                    df = pd.DataFrame(table_out,
                                      columns=['Gene',
                                               'Description'])
                    df.sort('Gene')
                    with open(output[0], "a") as out:
                        out.write("============================\n")
                        out.write("Target-site: " + str(prefix) + MOTIF + "\n")
                        out.write("Number-genes: " + str(num_genes) + "\n")
                        df.to_csv(out, sep='\t', index=False, quotechar='', quoting=csv.QUOTE_NONE, escapechar='\\')
                    table_out = [] #clear list
        # SECTION | CALCULATING SUMMARY STATISTICS ##################################
            total_genes = sum(num_genes_per_motif.values())
        # SECTION | DISPLAY SUMMARY OF FOUND GENE STATS #############################
            with open(output[1], "w") as out:
                out.write("Target-site\tNumber-genes\tPercentage\n")
                for k, v in num_genes_per_motif.items():
                    percent = round(v / total_genes, 4)
                    out.write(k + MOTIF + "\t" + str(v) + "\t" + str(percent) + "\n")

rule get_seqs_with_motif:
    input:
        "genome/READY.{A}.fa"
    params:
        motif = MOTIF
    output:
        "genome/{A}.motif.fa"
    shell:
        """
        for prefix in $(eval echo {{A,C,T,G}}{{A,C,T,G}}); do
            echo "*****BEGIN ${{prefix}}" >> {output}
            grep -i -B 1 "${{prefix}}{params}" {input} | sed '/^--$/d' >> {output}
            echo "*****END" >> {output}
        done
        """
