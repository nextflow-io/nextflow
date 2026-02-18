#!/bin/bash

exit_status=0

NXF_CMD=${NXF_CMD:-nextflow}
NXF_FILES=${*:-'*.nf'}

for pipeline in $NXF_FILES ; do

    echo "> Running test: $pipeline"

    # run the pipeline
    $NXF_CMD -q run "$pipeline" > .out

    # verify output (if ground truth exists)
    outfile="$(basename "$pipeline" .nf).out"

    if [[ -f "$outfile" ]] ; then
        sort "$outfile" > a.out
        sort .out > b.out
        diff a.out b.out
        status=$? ; [ $status -eq 0 ] || exit_status=$status
    else
        echo "> $outfile not found, skipping"
    fi

done

rm -rf .nextflow* work .out a.out b.out

[ $exit_status -eq 0 ] || false
