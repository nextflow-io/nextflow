/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

import java.lang.management.ManagementFactory

import com.sun.management.OperatingSystemMXBean
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.daemon.IgGridFactory
import nextflow.executor.IgBaseTask
import nextflow.file.FileHelper
import nextflow.util.MemoryUnit
import org.apache.ignite.Ignition
import org.apache.ignite.spi.collision.CollisionContext
import org.apache.ignite.spi.collision.CollisionJobContext

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class ResourceContext implements Serializable {

    static private transient OperatingSystemMXBean _bean

    private String _role

    private String _hostName

    private int _availCpus

    private MemoryUnit _availMemory

    private transient Collection<CollisionJobContext> activeJobs

    private transient Collection<CollisionJobContext> waitingJobs

    private int activeCount

    private int waitingCount

    private MemoryUnit freeMemory

    private int freeCpus


    @PackageScope
    ResourceContext() {

    }

    ResourceContext(CollisionContext ctx) {
        this(ctx.activeJobs(), ctx.waitingJobs());
    }

    ResourceContext(Collection<CollisionJobContext> activeJobs, Collection<CollisionJobContext> waitingJobs) {

        this.activeJobs = activeJobs;
        this.activeCount = activeJobs.size()

        this.waitingJobs = waitingJobs;
        this.waitingCount = waitingJobs.size()

        this.freeMemory = new MemoryUnit(getSystemMXBean().freePhysicalMemorySize)
        this.freeCpus = getAvailCpus()

        activeJobs.each { jobCtx ->

            if( jobCtx.job instanceof IgBaseTask ) {
                // count the number of used cpus
                final task = (IgBaseTask)jobCtx.job
                freeCpus -= task.resources.cpus
            }

        }

        // find out the actual free cpus
        if( freeCpus < 0 )
            freeCpus = 0

    }

    private static OperatingSystemMXBean getSystemMXBean() {
        if( !_bean ) {
            _bean = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
        }
        return _bean
    }

    /**
     * @return The node actual hostname
     */
    String getHostName() {
        if( _hostName ) {
            return _hostName
        }

        _hostName = System.getenv('HOSTNAME') ?: 'localhost'
    }


    /**
     * @return The total memory available
     */
    MemoryUnit getAvailMemory() {
        if( _availMemory ) {
            return _availMemory
        }

        _availMemory = new MemoryUnit(getSystemMXBean().getTotalPhysicalMemorySize())
    }

    /**
     * @return The {@link IgGridFactory#NODE_ROLE} attribute, either MASTER or WORKER
     */
    private String getRole() {
        if( _role ) {
            return _role
        }

        _role = Ignition.ignite(IgGridFactory.GRID_NAME).cluster().localNode().attribute(IgGridFactory.NODE_ROLE)
        log.trace "Local node role `$_role`"
        return _role
    }


    /**
     * @return The number of CPUs available
     */
    int getAvailCpus() {
        if( _availCpus ) {
            return _availCpus
        }

        _availCpus = getSystemMXBean().getAvailableProcessors()
        if( getRole() == IgGridFactory.ROLE_MASTER ) {
            // reserve some CPUs for nextflow/ignite scheduling activity on the master node
            if( _availCpus > 8 ) {
                _availCpus -= 2
            }
            else if( _availCpus > 3 ) {
                _availCpus -= 1
            }
        }
        return _availCpus
    }

    int getBusyCpus() {
        getAvailCpus() - freeCpus
    }

    int getFreeCpus() {
        freeCpus
    }

    MemoryUnit getFreeMemory() {
        freeMemory
    }

    /**
     * @return The actual free space in the node local storage
     */
    MemoryUnit getFreeDisk() {
        final free = FileHelper.getLocalTempPath().toFile().getFreeSpace()
        new MemoryUnit(free)
    }

    Collection<CollisionJobContext> getActiveJobs() { activeJobs }

    Collection<CollisionJobContext> getWaitingJobs() { waitingJobs }

    boolean canActivate( CollisionJobContext jobCtx )  {
        if( !(jobCtx.job instanceof IgBaseTask) ) {
            return true
        }

        final task = (IgBaseTask)jobCtx.job

        if( task.resources.disk > freeDisk ) {
            log.trace "Task waiting for disk storage > $task -- request: ${task.resources.disk}; free: ${freeDisk}"
            return false
        }

        if( task.resources.cpus > freeCpus ) {
            log.trace "Task waiting for cpus > $task -- request: ${task.resources.cpus}; free: ${freeCpus}"
            return false
        }

        if( task.resources.memory > freeMemory ) {
            log.trace "Task waiting for memory > $task -- request: ${task.resources.memory}; free: ${freeMemory}"
            return false
        }

        return true
    }

    void consumeTaskResource( CollisionJobContext jobCtx ) {
        final task = (IgBaseTask)jobCtx.job

        freeCpus -= task.resources.cpus
        freeMemory -= task.resources.memory
        if( log.isTraceEnabled() )
            log.trace "Task activated > $task -- Pending; ${waitingJobs.size()} (was: ${waitingCount}) - active: ${activeJobs.size()} (was: $activeCount)"

    }

    @Override
    String toString() {
        "cpus > avail: ${availCpus} free: ${freeCpus} - mem > avail: ${availMemory} free: $freeMemory} - disk > free: ${freeDisk}"
    }

}
