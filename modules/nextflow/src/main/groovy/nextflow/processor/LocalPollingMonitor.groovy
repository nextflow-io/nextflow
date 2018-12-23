/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.processor
import java.lang.management.ManagementFactory

import com.sun.management.OperatingSystemMXBean
import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.Duration
import nextflow.util.MemoryUnit

/**
 * Task polling monitor specialized for local execution. It manages tasks scheduling
 * taking into account task resources requests (cpus and memory)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class LocalPollingMonitor extends TaskPollingMonitor {

    static private OperatingSystemMXBean OS = { (OperatingSystemMXBean) ManagementFactory.getOperatingSystemMXBean() }()

    /**
     * Number of `free` CPUs available to execute pending tasks
     */
    private int availCpus

    /**
     * Total number of CPUs available in the system
     */
    private final int maxCpus

    /**
     * Amount of `free` memory available to execute pending tasks
     */
    private long availMemory

    /**
     * Total amount of memory available in the system
     */
    private final long maxMemory

    /**
     * Create the task polling monitor with the provided named parameters object.
     * <p>
     * Valid parameters are:
     * <li>name: The name of the executor for which the polling monitor is created
     * <li>session: The current {@code Session}
     * <li>capacity: The maximum number of this monitoring queue
     * <li>pollInterval: Determines how often a poll occurs to check for a process termination
     * <li>dumpInterval: Determines how often the executor status is written in the application log file
     *
     * @param params
     */
    protected LocalPollingMonitor(Map params) {
        super(params)
        this.availCpus = maxCpus = params.cpus as int
        this.availMemory = maxMemory = params.memory as long
        assert availCpus>0, "Local avail `cpus` attribute cannot be zero"
        assert availMemory>0, "Local avail `memory` attribute cannot zero"
    }

    /**
     * Creates an instance of {@link LocalPollingMonitor}
     *
     * @param session
     *      The current {@link Session} object
     * @param name
     *      The name of the executor that created this tasks monitor
     * @return
     *      An instance of {@link LocalPollingMonitor}
     */
    static LocalPollingMonitor create(Session session, String name) {
        assert session
        assert name

        final defPollInterval = Duration.of('100ms')
        final pollInterval = session.getPollInterval(name, defPollInterval)
        final dumpInterval = session.getMonitorDumpInterval(name)

        final int cpus = configCpus(session,name)
        final long memory = configMem(session,name)
        final int size = session.getQueueSize(name, OS.getAvailableProcessors())

        log.debug "Creating local task monitor for executor '$name' > cpus=$cpus; memory=${new MemoryUnit(memory)}; capacity=$size; pollInterval=$pollInterval; dumpInterval=$dumpInterval"

        new LocalPollingMonitor(
                name: name,
                cpus: cpus,
                memory: memory,
                session: session,
                capacity: size,
                pollInterval: pollInterval,
                dumpInterval: dumpInterval,
        )
    }

    @PackageScope
    static int configCpus(Session session, String name) {
        int cpus = session.getExecConfigProp(name, 'cpus', 0) as int

        if( !cpus )
            cpus = OS.getAvailableProcessors()

        return cpus
    }

    @PackageScope
    static long configMem(Session session, String name) {
        (session.getExecConfigProp(name, 'memory', OS.getTotalPhysicalMemorySize()) as MemoryUnit).toBytes()
    }

    /**
     * @param handler
     *      A {@link TaskHandler} instance
     * @return
     *      The number of cpus requested to execute the specified task
     */
    private static int cpus(TaskHandler handler) {
        handler.task.getConfig()?.getCpus()
    }

    /**
     *
     * @param handler
     *      A {@link TaskHandler} instance
     * @return
     *      The amount of memory (bytes) requested to execute the specified task
     */
    private static long mem(TaskHandler handler) {
        handler.task.getConfig()?.getMemory()?.toBytes() ?: 1
    }

    /**
     * Determines if a task can be submitted for execution checking if the resources required
     * (cpus and memory) match the amount of avail resource
     *
     * @param handler
     *      The {@link TaskHandler} representing the task to be executed
     * @return
     *      {@code true} if enough resources are available, {@code false} otherwise
     * @throws
     *      ProcessUnrecoverableException When the resource request exceed the total
     *      amount of resources provided by the underlying system e.g. task requires 10 cpus
     *      and the system provide 8 cpus
     *
     */
    @Override
    protected boolean canSubmit(TaskHandler handler) {

        final taskCpus = cpus(handler)
        if( taskCpus>maxCpus )
            throw new ProcessUnrecoverableException("Process requirement exceed available CPUs -- req: $taskCpus; avail: $maxCpus")

        final taskMemory = mem(handler)
        if( taskMemory>maxMemory)
            throw new ProcessUnrecoverableException("Process requirement exceed available memory -- req: ${new MemoryUnit(taskMemory)}; avail: ${new MemoryUnit(maxMemory)}")

        final result = super.canSubmit(handler) && taskCpus <= availCpus && taskMemory <= availMemory
        if( !result && log.isTraceEnabled( ) ) {
            log.trace "Task `${handler.task.name}` cannot be scheduled -- taskCpus: $taskCpus <= availCpus: $availCpus && taskMemory: ${new MemoryUnit(taskMemory)} <= availMemory: ${new MemoryUnit(availMemory)}"
        }
        return result
    }

    /**
     * Submits a task for execution allocating the resources (cpus and memory)
     * requested by the task
     * @param handler
     *      The {@link TaskHandler} representing the task to be executed
     */
    @Override
    protected void submit(TaskHandler handler) {
        super.submit(handler)
        availCpus -= cpus(handler)
        availMemory -= mem(handler)
    }

    /**
     * When a task completes its execution remove it from tasks polling queue
     * restoring the allocated resources i.e. cpus and memory
     *
     * @param handler
     *      The {@link TaskHandler} instance representing the task that completed its execution
     * @return
     *      {@code true} when the task is successfully removed from polling queue,
     *      {@code false} otherwise
     */
    protected boolean remove(TaskHandler handler) {
        final result = super.remove(handler)
        if( result ) {
            availCpus += cpus(handler)
            availMemory += mem(handler)
        }
        return result
    }
}
