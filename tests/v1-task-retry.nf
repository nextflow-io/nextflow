#!/usr/bin/env nextflow
// parser_v1: implicit environment variables

workflow {
  foo()
}

process foo {
    debug true
    time { 1.h * task.attempt }
    memory { 1.GB * task.attempt }
    errorStrategy { task.exitStatus == 5 && task.attempt<3 ? 'retry' : 'terminate' }
    maxErrors 10
    maxRetries 10

    script:
    """
    if [[ -f $PWD/marker ]]; then
    	echo DONE - mem: $task.memory - time: $task.time
    	exit 0
    else
    	echo FAIL
    	touch $PWD/marker
    	exit 5;
    fi
    """

}
