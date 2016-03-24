/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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


package nextflow.scheduler

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.MemoryUnit

/**
 * Model the computing resources required by a task
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString
@EqualsAndHashCode
@CompileStatic
class JobComputeResources implements Serializable {

    int cpus

    MemoryUnit memory

    MemoryUnit disk

    Duration time

    JobComputeResources() {}

    JobComputeResources( TaskRun task ) {
        cpus = task.config.getCpus()
        memory = task.config.getMemory()
        disk = task.config.getDisk()
        time = task.config.getTime()
    }

    @Override
    String toString() {
        "${this.class.simpleName} > cpus: ${cpus} - mem: ${memory} - disk: ${disk} - time: ${time}"
    }
}
