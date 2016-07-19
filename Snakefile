#!/usr/bin/env python

###############################################################################


###############################################################################

rule all:
    input:
        'genome/CDS.no_breaks.fa'

rule remove_fasta_linebreaks:
    input:
        'genome/Homo_sapiens.GRCh38.cds.all.fa'
    output:
        'genome/CDS.no_breaks.fa'
    shell:
        """
        awk '!/^>/ {{ printf "%s", $0; n = "\\n" }}
        /^>/ {{ print n $0; n = "" }}
        END {{ printf "%s", n }}
        ' {input} > {output}
        """
