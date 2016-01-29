package nextflow.collision
import java.lang.management.ManagementFactory

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.executor.IgBaseTask
import nextflow.util.MemoryUnit
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

    private String hostName

    private int availCpus

    private MemoryUnit availMemory

    private JobStealingCollisionSpi delegate

    CustomStealingCollisionSpi() {
        delegate = new JobStealingCollisionSpi()
        delegate.setActiveJobsThreshold(0)
        delegate.setWaitJobsThreshold(5)
        def field = JobStealingCollisionSpi.getDeclaredField('log')
        field.setAccessible(true)
        field.set(delegate, new Slf4jLogger().getLogger(JobStealingCollisionSpi))
    }

    @Override
    String getHostName() {
        return hostName
    }

    @Override
    int getAvailCpus() {
        return availCpus
    }

    @Override
    MemoryUnit getAvailMemory() {
        return availMemory
    }

    @Override
    void onCollision(CollisionContext ctx) {

        final bean = (com.sun.management.OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
        final freeMemory = new MemoryUnit(bean.freePhysicalMemorySize)

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

        final bean = (com.sun.management.OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean()
        this.availCpus = bean.availableProcessors
        this.availMemory = new MemoryUnit(bean.totalPhysicalMemorySize)
        this.hostName = System.getenv('HOSTNAME') ?: 'localhost'

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
