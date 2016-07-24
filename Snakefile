#!/usr/bin/env python

"""
Author: K. Chen
Affiliation: NCI/NIH
Aim: A pipeline to find and analzye possible non-canonical miR targets using simple seed matching.
Date: Mon Nov  2 14:03:11 CET 2015
Run: snakemake
"""

import re
import pandas as pd

###############################################################################

NAMES = ['cds', 'cdna']
MOTIF = "GCGGAAC"

###############################################################################

rule all:
    input:
        expand("results/{A}.txt", A=NAMES)

rule get_analysis_of_seqs:
    input:
        "genome/{A}.no_breaks.canon.motif.fa"
    output:
        "results/{A}.txt"
    run:
        with open(str(input), "rt") as input_file:
            table_out = []
            num_genes = 0
            prefix = ""
            for line in input_file:
                if "*****BEGIN" in line:
                    num_genes = 0  # reset
                    prefix = line.split(" ")[1].rstrip()

                # Collate information
                if "gene_symbol:" in line:
                    num_genes += 1
                    try:
                        gene_symbol = line.split("gene_symbol:")[1].split(" description:")[0]
                    except IndexError:
                        gene_symbol = "null"
                    try:
                        gene_desc = line.split("description:")[1].split(" [")[0]
                    except IndexError:
                        gene_desc = "null"
                    table_out.append([gene_symbol, gene_desc])

                if "*****END" in line:
                    df = pd.DataFrame(table_out,
                                      columns=["Gene",
                                               "Description"])
                    with open(output[0], "a") as out:
                        out.write("============================\n")
                        out.write("Target-site: " + str(prefix) + MOTIF + "\n")
                        out.write("Number-genes: " + str(num_genes) + "\n\n")
                        df.to_csv(out, sep='\t', index=False)

rule get_seqs_with_motif:
    input:
        "genome/{A}.no_breaks.canon.fa"
    params:
        motif = MOTIF
    output:
        "genome/{A}.no_breaks.canon.motif.fa"
    shell:
        """
        for prefix in $(eval echo {{A,C,T,G}}{{A,C,T,G}}); do
            echo "*****BEGIN ${{prefix}}" >> {output}
            grep -i -B 1 "${{prefix}}{params}" {input} >> {output}
            echo "*****END" >> {output}
        done
        """

rule get_longest_seq:
    input:
        "genome/{A}.no_breaks.fa"
    output:
        "genome/{A}.no_breaks.canon.fa"
    shell:
        """
        python longest_seq.py {input} > {output}
        """

rule remove_fasta_linebreaks:
    input:
        "genome/Homo_sapiens.GRCh38.{a}.all.fa"
    output:
        "genome/{a}.no_breaks.fa"
    shell:
        """
        awk '!/^>/ {{ printf "%s", $0; n = "\\n" }}
        /^>/ {{ print n $0; n = "" }}
        END {{ printf "%s", n }}
        ' {input} > {output}
        """

rule download_cdna:
    output:
        "genome/Homo_sapiens.GRCh38.cdna.all.fa"
    shell:
        "curl 'ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/cdna//Homo_sapiens.GRCh38.cdna.all.fa.gz' | gunzip -c > {output}"

rule download_cds:
    output:
        "genome/Homo_sapiens.GRCh38.cds.all.fa"
    shell:
        "curl 'ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/cds//Homo_sapiens.GRCh38.cds.all.fa.gz' | gunzip -c > {output}"
