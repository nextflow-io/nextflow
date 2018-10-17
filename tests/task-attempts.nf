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


process foo {

    memory { def x=task.attempt * 1.GB; println "attempt: $task.attempt; memory: $x"; return x }
    time { def x=1.h * task.attempt; println "attempt: $task.attempt; time: $x"; return x }
    errorStrategy { task.exitStatus == 5 && task.attempt<3 ? 'retry' : 'terminate' }
    maxErrors 10
    maxRetries 10

    script:
    """
    exit 5
    """

}
