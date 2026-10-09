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

import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicLong

import groovy.transform.CompileStatic
import nextflow.exception.ProcessUnrecoverableException
import nextflow.util.MemoryUnit

/**
 * Keeps track of the cpus and memory requested by running tasks,
 * relative to a fixed total. A total of zero means no limit for
 * that resource.
 */
@CompileStatic
class ResourceTracker {

    final int totalCpus

    final long totalMemory

    private final AtomicInteger availCpus

    private final AtomicLong availMemory

    ResourceTracker(int totalCpus, long totalMemory) {
        assert totalCpus >= 0, "Executor `cpus` setting cannot be a negative value"
        assert totalMemory >= 0, "Executor `memory` setting cannot be a negative value"
        this.totalCpus = totalCpus
        this.totalMemory = totalMemory
        this.availCpus = new AtomicInteger(totalCpus)
        this.availMemory = new AtomicLong(totalMemory)
    }

    int availableCpus() { availCpus.get() }

    long availableMemory() { availMemory.get() }

    /**
     * Get the number of cpus requested by a task, or the sum
     * of its children for a job array.
     *
     * @param handler
     */
    static int cpus(TaskHandler handler) {
        if( handler.task instanceof TaskArrayRun ) {
            int result = 0
            for( TaskHandler child : ((TaskArrayRun)handler.task).children )
                result += cpus(child)
            return result
        }
        handler.task.getConfig()?.getCpus() ?: 1
    }

    /**
     * Get the amount of memory (bytes) requested by a task, or the sum
     * of its children for a job array.
     *
     * @param handler
     */
    static long memory(TaskHandler handler) {
        if( handler.task instanceof TaskArrayRun ) {
            long result = 0
            for( TaskHandler child : ((TaskArrayRun)handler.task).children )
                result += memory(child)
            return result
        }
        handler.task.getConfig()?.getMemory()?.toBytes() ?: 1L
    }

    /**
     * Fail if a task requests more resources than the total.
     *
     * @param handler
     */
    void validate(TaskHandler handler) {
        final array = handler.task instanceof TaskArrayRun
        final prefix = array ? 'Array' : 'Task'
        final suffix = array ? " (array size: ${((TaskArrayRun)handler.task).children.size()})" : ''

        final taskCpus = cpus(handler)
        if( totalCpus && taskCpus > totalCpus )
            throw new ProcessUnrecoverableException("$prefix requirement exceeds available CPUs -- req: $taskCpus$suffix; avail: $totalCpus")

        final taskMemory = memory(handler)
        if( totalMemory && taskMemory > totalMemory )
            throw new ProcessUnrecoverableException("$prefix requirement exceeds available memory -- req: ${new MemoryUnit(taskMemory)}$suffix; avail: ${new MemoryUnit(totalMemory)}")
    }

    boolean canAcquire(TaskHandler handler) {
        (!totalCpus || cpus(handler) <= availCpus.get()) && (!totalMemory || memory(handler) <= availMemory.get())
    }

    void acquire(TaskHandler handler) {
        if( totalCpus )
            availCpus.addAndGet(-cpus(handler))
        if( totalMemory )
            availMemory.addAndGet(-memory(handler))
    }

    /**
     * Release the resources of a task. Each child of a job array
     * releases its own share of the array's resources.
     *
     * @param handler
     */
    void release(TaskHandler handler) {
        if( totalCpus )
            availCpus.addAndGet(cpus(handler))
        if( totalMemory )
            availMemory.addAndGet(memory(handler))
    }

    String toString() {
        "ResourceTracker[cpus=${availCpus.get()}/$totalCpus; memory=${new MemoryUnit(availMemory.get())}/${new MemoryUnit(totalMemory)}]"
    }
}
