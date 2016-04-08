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

package nextflow.daemon

import groovy.transform.CompileStatic
import nextflow.util.MemoryUnit
import org.apache.ignite.mxbean.MXBeanDescription
import org.apache.ignite.spi.IgniteSpiManagementMBean

/**
 * MBean interface required to expose SPI properties in the Ignite config console
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
interface CustomStealingCollisionSpiMBean extends IgniteSpiManagementMBean {

    @MXBeanDescription("Compute host name")
    String getHostName()

    @MXBeanDescription("Total cpus available")
    int getAvailCpus()

    @MXBeanDescription("Total memory available")
    MemoryUnit getAvailMemory()

}