#!/usr/bin/env python

###############################################################################


###############################################################################

rule all:
    input:
        expand('results/tabular/{A}.isomir.tsv', A=SAMPLES),
        expand('results/text/{A}.mirna.txt', A=SAMPLES)

rule get:
    input:
        'Homo_sapiens.GRCh38.cds.all.fa'
    output:
        'CDS_no_linebreaks.fa'
    shell:
        """
        awk '/^>/{print s? s'\n'$0:$0;s='';next}{s=s sprintf('%s',$0)}END{if(s)print s}' > {output}
        """
