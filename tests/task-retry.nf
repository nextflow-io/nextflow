#!/usr/bin/env nextflow
/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

echo true

process foo {
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
