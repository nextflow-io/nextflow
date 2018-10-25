#!/bin/bash -ue
# define box inputs
CONT_INPUT_FASTA=sample.fa
CONT_OUTPUT_FILE=result.fa
# launch box run
docker run -i -e "CONT_INPUT_FASTA=sample.fa" -e "CONT_OUTPUT_FILE=result.fa" -v /Users/olafurh/git/nextflow/nextflow-private/tests:/Users/olafurh/git/nextflow/nextflow-private/tests -v "$PWD":"$PWD" -w "$PWD" --name $NXF_BOXID nextflow/tcoffee
