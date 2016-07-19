#!/usr/bin/env python

"""
Author: K. Chen
Affiliation: NCI/NIH
Aim: A pipeline to find and analzye possible non-canonical miR targets using simple seed matching.
Date: Mon Nov  2 14:03:11 CET 2015
Run: snakemake
"""
###############################################################################

###############################################################################

rule all:
    input:
        'genome/CDS.no_breaks.fa',
        'genome/CDNA.no_breaks.fa'

rule find_motif:
    shell:
        # Find all sequences in sample that match motif (w 2nt permutation
        # prefix
        for prefix in $(eval echo {A, C, T, G}{A, C, T, G}); do
            grep - i - B 1 "${prefix}GCGGAAC" CDS_no_linebreaks.fa > CDS.${prefix}GCGGAAC.fa
        done

        # Find all unique genes in matched motifs
        for prefix in $(eval echo {A, C, T, G}{A, C, T, G}); do
            grep - o 'gene_symbol:\S*' CDS.${prefix}GCGGAAC.fa | sort | uniq - u > GENES.${prefix}GCGGAAC.fa
        done


# When alternative transcript isoforms are reported for same Ensembl gene
# ID, only longest coding sequence retained
rule get_longest_cds:

rule get_cdna:
    output:
        'genome/Homo_sapiens.GRCh38.cdna.all.fa'
    shell:
        "curl 'ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/cdna//Homo_sapiens.GRCh38.cdna.all.fa.gz' | gunzip -c > {output}"

rule get_coding_sequences:
    output:
        'genome/Homo_sapiens.GRCh38.cds.all.fa'
    shell:
        "curl 'ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/cds//Homo_sapiens.GRCh38.cds.all.fa.gz' | gunzip -c > {output}"

rule remove_fasta_linebreaks:
    input:
        'genome/Homo_sapiens.GRCh38.cds.all.fa',
        'genome/Homo_sapiens.GRCh38.cdna.all.fa'
    output:
        'genome/CDS.no_breaks.fa',
        'genome/CDNA.no_breaks.fa'
    shell:
        """
        awk '!/^>/ {{ printf "%s", $0; n = "\\n" }}
        /^>/ {{ print n $0; n = "" }}
        END {{ printf "%s", n }}
        ' {input[0]} > {output[0]}

        awk '!/^>/ {{ printf "%s", $0; n = "\\n" }}
        /^>/ {{ print n $0; n = "" }}
        END {{ printf "%s", n }}
        ' {input[1]} > {output[1]}
        """
