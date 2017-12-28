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

package nextflow.cloud.types

import groovy.transform.Immutable
import nextflow.util.MemoryUnit

/**
 * Models a cloud instance type
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Immutable
class CloudInstanceType implements Serializable, Cloneable {

    /**
     * Instance type ID as assigned by the cloud provider
     */
    String id

    /**
     * Number of CPUs
     */
    int cpus

    /**
     * Amount of memory (RAM)
     */
    MemoryUnit memory

    /**
     * Amount of local storage
     */
    MemoryUnit disk

    /**
     * Number of disks
     */
    int numOfDisks

    String toString() {
        "id=$id; cpus=$cpus; mem=$memory; disk=$disk x $numOfDisks"
    }
}
