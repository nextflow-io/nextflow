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
import java.lang.management.ManagementFactory

import com.sun.management.OperatingSystemMXBean
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Global
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

    private transient Collection<CollisionJobContext> activeJobs

    private transient Collection<CollisionJobContext> waitingJobs

    private transient String message

    private int activeCount

    private int waitingCount

    private MemoryUnit freeMemory

    private MemoryUnit freeDisk

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

        this.freeDisk = getAvailDisk()
        this.freeMemory = getAvailMemory()
        this.freeCpus = getAvailCpus()

        activeJobs.each { jobCtx ->

            if( jobCtx.job instanceof IgBaseTask ) {
                // count the number of used cpus
                final task = (IgBaseTask)jobCtx.job
                freeCpus -= task.resources.cpus
                freeMemory -=  task.resources.memory
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
    @Memoized
    String getHostName() {
        System.getenv('HOSTNAME') ?: 'localhost'
    }


    /**
     * @return The total memory available
     */
    @Memoized
    MemoryUnit getAvailMemory() {
        new MemoryUnit(getSystemMXBean().getTotalPhysicalMemorySize())
    }

    /**
     * @return The {@link IgGridFactory#NODE_ROLE} attribute, either MASTER or WORKER
     */
    @Memoized
    String getRole() {
        Ignition.ignite(IgGridFactory.GRID_NAME).cluster().localNode().attribute(IgGridFactory.NODE_ROLE)
    }


    /**
     * @return The number of CPUs available
     */
    @Memoized
    int getAvailCpus() {

        final int result = getSystemMXBean().getAvailableProcessors()
        if( getRole() == IgGridFactory.ROLE_MASTER ) {
            // reserve some CPUs for nextflow/ignite scheduling activity on the master node
            result -= getMasterMinCpus()
        }
        log.trace1 "Avail cpus: $result"
        return result
    }

    @CompileDynamic
    private int getMasterMinCpus() {
        Global.session?.config?.cluster?.masterNodeMinCpus ?: 0
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

    MemoryUnit getFreeDisk() {
        freeDisk
    }

    /**
     * @return The actual free space in the node local storage
     */
    MemoryUnit getAvailDisk() {
        final free = FileHelper.getLocalTempPath().toFile().getFreeSpace()
        new MemoryUnit(free)
    }

    Collection<CollisionJobContext> getActiveJobs() { activeJobs }

    Collection<CollisionJobContext> getWaitingJobs() { waitingJobs }

    String getMessage() {
        message
    }

    /**
     * Determine if a task can be activated depending the computing resources requested and
     * the ones available.
     *
     * @param jobCtx
     * @return
     */
    boolean canActivate( CollisionJobContext jobCtx )  {
        if( jobCtx.job instanceof IgBaseTask ) {
            final task = (IgBaseTask)jobCtx.job

            if( task.resources.cpus > freeCpus ) {
                this.message = "task `$task.taskId` exceed avail cpus > requested: ${task.resources.cpus}; free: ${freeCpus}"
                return false
            }

            if( task.resources.memory > freeMemory ) {
                this.message = "task `$task.taskId` exceed avail memory > requested: ${task.resources.memory}; free: ${freeMemory}"
                return false
            }

            if( task.resources.disk > freeDisk ) {
                this.message = "task `$task.taskId` exceed avail storage > requested: ${task.resources.disk}; free: ${freeDisk}"
                return false
            }
        }

        return true
    }

    /**
     * Consume the job resources
     *
     * @param jobCtx
     */
    void consumeTaskResource( CollisionJobContext jobCtx ) {

        if( jobCtx.getJob() instanceof IgBaseTask ) {
            final task = (IgBaseTask) jobCtx.job
            freeCpus -= task.resources.cpus
            freeMemory -= task.resources.memory
            freeDisk -= task.resources.disk
        }
        else {
            freeCpus -= 1
        }

        log.trace1("Task activated > ${jobCtx.job} -- Jobs pending: ${waitingJobs.size()} (was: ${waitingCount}) - active: ${activeJobs.size()} (was: $activeCount)")
    }

    @Override
    String toString() {
        "cpus > tot: ${availCpus} free: ${freeCpus} - memory > tot: ${availMemory} free: $freeMemory} - disk > free: ${freeDisk}"
    }

}
