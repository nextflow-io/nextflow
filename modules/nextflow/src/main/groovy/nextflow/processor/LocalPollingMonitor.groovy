/*
 * Copyright 2013-2026, Seqera Labs
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
import nextflow.executor.ExecutorConfig
import nextflow.exception.ProcessUnrecoverableException
import nextflow.executor.local.LocalTaskHandler
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
     * Tracks the total and available accelerators in the system
     */
    private AcceleratorTracker acceleratorTracker

    /**
     * Create the task polling monitor with the provided named parameters object.
     * <p>
     * Valid parameters are:
     * <li>name: The name of the executor for which the polling monitor is created
     * <li>session: The current {@code Session}
     * <li>config: The `executor` configuration settings
     * <li>capacity: The maximum number of this monitoring queue
     * <li>pollInterval: Determines how often a poll occurs to check for a process termination
     * <li>dumpInterval: Determines how often the executor status is written in the application log file
     *
     * @param params
     */
    protected LocalPollingMonitor(Map params) {
        super(params)
        final cpus = params.cpus as int
        final memory = params.memory as long
        assert cpus>0, "Local avail `cpus` attribute cannot be zero"
        assert memory>0, "Local avail `memory` attribute cannot zero"
        this.resourceTracker = new ResourceTracker(cpus, memory)
        this.acceleratorTracker = AcceleratorTracker.create()
    }

    /**
     * Creates an instance of {@link LocalPollingMonitor}
     *
     * @param session
     *      The current {@link Session} object
     * @param config
     *      The `executor` configuration settings
     * @param name
     *      The name of the executor that created this tasks monitor
     * @return
     *      An instance of {@link LocalPollingMonitor}
     */
    static LocalPollingMonitor create(Session session, ExecutorConfig config, String name) {
        assert session
        assert config
        assert name

        final pollInterval = config.getPollInterval(name, Duration.of('100ms'))
        final dumpInterval = config.getMonitorDumpInterval(name)
        final cpus = configCpus(config, name)
        final memory = configMem(config, name)
        final size = config.getQueueSize(name, OS.getAvailableProcessors())

        log.debug "Creating local task monitor for executor '$name' > cpus=$cpus; memory=${new MemoryUnit(memory)}; capacity=$size; pollInterval=$pollInterval; dumpInterval=$dumpInterval"

        new LocalPollingMonitor(
                name: name,
                cpus: cpus,
                memory: memory,
                session: session,
                config: config,
                capacity: size,
                pollInterval: pollInterval,
                dumpInterval: dumpInterval,
        )
    }

    @PackageScope
    static int configCpus(ExecutorConfig config, String name) {
        int cpus = config.getExecConfigProp(name, 'cpus', 0) as int

        if( !cpus )
            cpus = OS.getAvailableProcessors()

        return cpus
    }

    @PackageScope
    static long configMem(ExecutorConfig config, String name) {
        final memory = config.getExecConfigProp(name, 'memory', OS.getTotalPhysicalMemorySize()) as MemoryUnit
        return memory.toBytes()
    }

    /**
     * @param handler
     *      A {@link TaskHandler} instance
     * @return
     *      The number of accelerators requested to execute the specified task
     */
    private static int accelerators(TaskHandler handler) {
        handler.task.getConfig()?.getAccelerator()?.getRequest() ?: 0
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
        final taskAccelerators = accelerators(handler)
        if( acceleratorTracker.name() != null && taskAccelerators > acceleratorTracker.total() )
            throw new ProcessUnrecoverableException("Process requirement exceeds available accelerators -- req: $taskAccelerators; avail: ${acceleratorTracker.total()}")

        final accelOk = acceleratorTracker.name() == null || taskAccelerators <= acceleratorTracker.available()
        final result = super.canSubmit(handler) && accelOk
        if( !result && log.isTraceEnabled( ) ) {
            log.trace "Task `${handler.task.name}` cannot be scheduled -- taskCpus: ${ResourceTracker.cpus(handler)} <= availCpus: ${resourceTracker.availableCpus()} && taskMemory: ${new MemoryUnit(ResourceTracker.memory(handler))} <= availMemory: ${new MemoryUnit(resourceTracker.availableMemory())} && taskAccelerators: $taskAccelerators <= availAccelerators: ${acceleratorTracker.name() != null ? acceleratorTracker.available() : 'n/a'}"
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
        final taskAccelerators = accelerators(handler)
        if( handler instanceof LocalTaskHandler && acceleratorTracker.name() != null && taskAccelerators > 0 ) {
            handler.acceleratorEnv = acceleratorTracker.name()
            handler.acceleratorIds = acceleratorTracker.acquire(taskAccelerators)
        }

        try {
            super.submit(handler)
        }
        catch( Throwable e ) {
            if( handler instanceof LocalTaskHandler && handler.acceleratorIds )
                acceleratorTracker.release(handler.acceleratorIds)
            throw e
        }
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
    @Override
    protected boolean remove(TaskHandler handler) {
        final result = super.remove(handler)
        if( result && handler instanceof LocalTaskHandler )
            acceleratorTracker.release(handler.acceleratorIds ?: Collections.<String>emptyList())
        return result
    }
}
