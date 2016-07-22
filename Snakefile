#!/usr/bin/env python

"""
Author: K. Chen
Affiliation: NCI/NIH
Aim: A pipeline to find and analzye possible non-canonical miR targets using simple seed matching.
Date: Mon Nov  2 14:03:11 CET 2015
Run: snakemake
"""

###############################################################################

NAMES = ['cds', 'cdna']

###############################################################################

rule all:
    input:
        expand('genome/{A}.no_breaks.canon.fa', A=NAMES)

rule get_longest_seq:
    input:
        'genome/{A}.no_breaks.fa'
    output:
        'genome/{A}.no_breaks.canon.fa'
    shell:
        """
        python longest_seq.py {input} > {output}
        """

rule remove_fasta_linebreaks:
    input:
        'genome/Homo_sapiens.GRCh38.{a}.all.fa'
    output:
        'genome/{a}.no_breaks.fa'
    shell:
        """
        awk '!/^>/ {{ printf "%s", $0; n = "\\n" }}
        /^>/ {{ print n $0; n = "" }}
        END {{ printf "%s", n }}
        ' {input} > {output}
        """

rule download_cdna:
    output:
        'genome/Homo_sapiens.GRCh38.cdna.all.fa'
    shell:
        "curl 'ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/cdna//Homo_sapiens.GRCh38.cdna.all.fa.gz' | gunzip -c > {output}"

rule download_cds:
    output:
        'genome/Homo_sapiens.GRCh38.cds.all.fa'
    shell:
        "curl 'ftp://ftp.ensembl.org/pub/release-85/fasta/homo_sapiens/cds//Homo_sapiens.GRCh38.cds.all.fa.gz' | gunzip -c > {output}"
