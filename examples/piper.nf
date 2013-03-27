#!/bin/env nextflow

echo true

query = '/Users/ptommaso/workspace/piper/tutorial/5_RNA_queries.fa'
dbpath = '/Users/ptommaso/workspace/piper/tutorial'
genomes = ['briggsae', 'elegans']


outFmt = '6 qseqid sseqid evalue score qgi bitscore length nident positive mismatch pident ppos qacc gaps gaopen qaccver qlen qframe qstart qend sframe sstart send'


blastResult = queue()

task ('blast') {
    input genomes
    output 'result_*': blastResult

    """
    blastn -db $dbpath/$genomes/blastdb/db -query $query -outfmt '$outFmt' > result_$genomes
    """

}



task ('exonerate') {
    input blastResult

    """
    ID=`basename $blastResult | sed 's/result_//'`
    exonerateRemapping.pl -query $query -mf2 $blastResult -targetGenomeFolder $dbpath/$ID/chr/ -exonerate_lines_mode 1000 -exonerate_success_mode 1 -ner no
    """

}