package nextflow.daemon
import java.lang.management.ManagementFactory

import com.sun.management.OperatingSystemMXBean
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.IgBaseTask
import nextflow.file.FileHelper
import nextflow.util.MemoryUnit
import org.apache.ignite.Ignition
import org.apache.ignite.logger.slf4j.Slf4jLogger
import org.apache.ignite.spi.IgniteSpiAdapter
import org.apache.ignite.spi.IgniteSpiException
import org.apache.ignite.spi.collision.CollisionContext
import org.apache.ignite.spi.collision.CollisionExternalListener
import org.apache.ignite.spi.collision.CollisionSpi
import org.apache.ignite.spi.collision.jobstealing.JobStealingCollisionSpi
import org.jetbrains.annotations.Nullable

/**
 * Extends stock {@link CustomStealingCollisionSpi} adding the ability to manage job requested resources
 * (cpus, memory, storage, etc)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CustomStealingCollisionSpi extends IgniteSpiAdapter implements CollisionSpi, CustomStealingCollisionSpiMBean {

    private String fHostName

    private int fAvailCpus

    private MemoryUnit fAvailMemory

    private JobStealingCollisionSpi delegate

    private String fRole

    private OperatingSystemMXBean fBean

    private boolean resourcesLogged

    CustomStealingCollisionSpi() {
        delegate = new JobStealingCollisionSpi()
        delegate.setActiveJobsThreshold(0)
        delegate.setWaitJobsThreshold(5)
        def field = JobStealingCollisionSpi.getDeclaredField('log')
        field.setAccessible(true)
        field.set(delegate, new Slf4jLogger().getLogger(JobStealingCollisionSpi))
    }

    private OperatingSystemMXBean getSystemMXBean() {
        if( !fBean ) {
            fBean = (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
        }
        return fBean
    }

    @Override
    String getHostName() {
        if( fHostName ) {
            return fHostName
        }

        fHostName = System.getenv('HOSTNAME') ?: 'localhost'
    }

    @Override
    int getAvailCpus() {
        if( fAvailCpus ) {
            return fAvailCpus
        }

        fAvailCpus = getSystemMXBean().getAvailableProcessors()
        if( getRole() == IgGridFactory.ROLE_MASTER ) {
            // reserve some CPUs for nextflow/ignite scheduling activity on the master node
            if( fAvailCpus > 8 ) {
                fAvailCpus -= 2
            }
            else if( fAvailCpus > 3 ) {
                fAvailCpus -= 1
            }
        }
        return fAvailCpus
    }

    @Override
    MemoryUnit getAvailMemory() {
        if( fAvailMemory ) {
            return fAvailMemory
        }

        fAvailMemory = new MemoryUnit(getSystemMXBean().getTotalPhysicalMemorySize())
    }

    private String getRole() {
        if( fRole ) {
            return fRole
        }

        fRole = Ignition.ignite(IgGridFactory.GRID_NAME).cluster().localNode().attribute(IgGridFactory.NODE_ROLE)
        log.trace "Local node role `$fRole`"
        return fRole
    }

    private logResources() {
        if( !resourcesLogged ) {
            resourcesLogged = true
            log.debug "Computing resources for node: `$hostName` [${getRole()}] > cpus: ${availCpus}; mem: ${availMemory}; free disk: ${freeDisk} (${FileHelper.getLocalTempPath()})"
        }
    }

    private MemoryUnit getFreeDisk() {
        final free = FileHelper.getLocalTempPath().toFile().getFreeSpace()
        new MemoryUnit(free)
    }

    @Override
    void onCollision(CollisionContext ctx) {

        logResources()

        final freeMemory = new MemoryUnit(getSystemMXBean().freePhysicalMemorySize)
        final activeJobs = ctx.activeJobs().size()
        final waitingJobs = ctx.waitingJobs().size()

        /*
         * get the count of used and free cpus
         */
        int busyCpus = 0
        ctx
            .activeJobs()
            .findAll{ it.job instanceof IgBaseTask }
            .collect { (IgBaseTask)it.job }
            .each { IgBaseTask job -> busyCpus += job.resources.cpus }

        final freeCpus = availCpus>busyCpus ? availCpus-busyCpus : 0

        if( log.isTraceEnabled() ) {
            log.trace "Node `$hostName` resources > cpus: $availCpus ($freeCpus) - mem: $availMemory ($freeMemory) - active: $activeJobs - waiting: $waitingJobs"
        }

        // activate waiting jobs that match avail resources
        ctx .waitingJobs() .each { jobCtx ->
            def task = (IgBaseTask)jobCtx.job

            if( task.resources.cpus > availCpus ) {
                if(log.isTraceEnabled()) log.trace "Task rejected - it requires more cpus than the available ones: ${task.resources.cpus} (${availCpus}) "
                jobCtx.cancel()
            }

            if( task.resources.memory > availMemory ) {
                if(log.isTraceEnabled()) log.trace "Task rejected - it requires more memory than the available one: ${task.resources.memory} (${availMemory}) "
                jobCtx.cancel()
            }

            if( task.resources.disk > freeDisk ) {
                if(log.isTraceEnabled()) log.trace "Task rejected - it requires more disk storage than the available one: ${task.resources.disk} (${freeDisk}) "
                jobCtx.cancel()
            }

            if( task.resources.cpus > freeCpus ) return
            if( task.resources.memory > freeMemory ) return

            if( jobCtx.activate() ) {
                freeCpus -= task.resources.cpus
                freeMemory -= task.resources.memory
                if( log.isTraceEnabled() )
                    log.trace "Activated task > $task -- pending ${ctx.waitingJobs().size()} (was ${waitingJobs}) - active ${ctx.activeJobs().size()} (was $activeJobs)"
            }
            else if(log.isTraceEnabled())  {
                log.trace "Failed to activate task > $task"
            }

        }

        // fallback on the job stealing strategy
        delegate.onCollision(ctx)
    }

    @Override
    void setExternalCollisionListener(@Nullable CollisionExternalListener listener) {
        // No-op.
    }

    /** {@inheritDoc} */
    @Override
    void spiStart(String gridName) throws IgniteSpiException {

        // Start SPI start stopwatch.
        startStopwatch();
        registerMBean(gridName, this, CustomStealingCollisionSpiMBean.class);

        log.debug(startInfo())
    }

    /** {@inheritDoc} */
    @Override
    void spiStop() throws IgniteSpiException {
        unregisterMBean();
        log.debug(stopInfo())
    }

}
