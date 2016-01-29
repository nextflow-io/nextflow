package nextflow.collision

import groovy.transform.CompileStatic
import nextflow.util.MemoryUnit
import org.apache.ignite.mxbean.MXBeanDescription
import org.apache.ignite.spi.IgniteSpiManagementMBean

/**
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