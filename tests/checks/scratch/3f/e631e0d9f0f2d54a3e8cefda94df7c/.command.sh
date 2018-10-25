#!/bin/bash -ue
if [[ -f /Users/olafurh/git/nextflow/nextflow-private/tests/checks/task-retry.nf/marker ]]; then
	echo DONE - mem: 2 GB - time: 2h
	exit 0
else
	echo FAIL
	touch /Users/olafurh/git/nextflow/nextflow-private/tests/checks/task-retry.nf/marker
	exit 5;
fi
