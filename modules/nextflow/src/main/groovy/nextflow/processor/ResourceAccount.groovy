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
import groovy.transform.PackageScope
import nextflow.util.MemoryUnit
/**
 * Tracks the amount of CPUs and memory reserved by the tasks submitted
 * by a task monitor, relative to a fixed maximum.
 *
 * The account is used to throttle the task submission when the `executor.cpus`
 * and/or `executor.memory` settings are specified.
 *
 */
@CompileStatic
class ResourceAccount {

    /**
     * Total number of CPUs allowed by the account, zero when unlimited
     */
    final int maxCpus

    /**
     * Total amount of memory in bytes allowed by the account, zero when unlimited
     */
    final long maxMemory

    /**
     * Number of CPUs not yet reserved
     */
    private final AtomicInteger availCpus

    /**
     * Amount of memory in bytes not yet reserved
     */
    private final AtomicLong availMemory

    /**
     * Create an account with unlimited resources
     */
    ResourceAccount() {
        this(0, 0)
    }

    ResourceAccount(int maxCpus, long maxMemory) {
        assert maxCpus >= 0, "Executor `cpus` setting cannot be a negative value"
        assert maxMemory >= 0, "Executor `memory` setting cannot be a negative value"
        this.maxCpus = maxCpus
        this.maxMemory = maxMemory
        this.availCpus = new AtomicInteger(maxCpus)
        this.availMemory = new AtomicLong(maxMemory)
    }

    boolean isUnlimited() { maxCpus == 0 && maxMemory == 0 }

    int availableCpus() { availCpus.get() }

    long availableMemory() { availMemory.get() }

    /**
     * @param handler A {@link TaskHandler} for the task
     * @return The number of CPUs requested by the specified task handler. Job arrays
     *      require the sum of the resources requested by all the array children.
     */
    static int cpusOf(TaskHandler handler) {
        if( handler.task instanceof TaskArrayRun ) {
            int result = 0
            for( TaskHandler child : ((TaskArrayRun)handler.task).children )
                result += child.task.config?.getCpus() ?: 1
            return result
        }
        handler.task.config?.getCpus() ?: 1
    }

    /**
     * @param handler A {@link TaskHandler} for the task
     * @return The amount of memory in bytes requested by the specified task handler.
     *      Job arrays require the sum of the resources requested by all the array children.
     */
    static long memOf(TaskHandler handler) {
        if( handler.task instanceof TaskArrayRun ) {
            long result = 0
            for( TaskHandler child : ((TaskArrayRun)handler.task).children )
                result += child.task.config?.getMemory()?.toBytes() ?: 0
            return result
        }
        handler.task.config?.getMemory()?.toBytes() ?: 0
    }

    /**
     * @return {@code true} if the resources requested by the specified task handler
     *      can be satisfied by the amount of resources remaining in this account
     */
    boolean canReserve(TaskHandler handler) {
        canReserve(cpusOf(handler), memOf(handler))
    }

    /**
     * @return {@code true} if the specified amount of resources can be satisfied
     *      by the amount of resources remaining in this account
     */
    boolean canReserve(int cpus, long memory) {
        if( isUnlimited() )
            return true
        if( maxCpus > 0 && (cpus > maxCpus || cpus > availCpus.get()) )
            return false
        if( maxMemory > 0 && (memory > maxMemory || memory > availMemory.get()) )
            return false
        return true
    }

    /**
     * Reserve the resources requested by the specified task handler.
     *
     * @return A reservation object to be passed to {@link #release(ResourceAccount.Reservation)}
     *      when the task is removed from the running queue, or {@code null} when the
     *      account is unlimited i.e. no tracking is required
     */
    Reservation reserve(TaskHandler handler) {
        if( isUnlimited() )
            return null
        final result = new Reservation(cpusOf(handler), memOf(handler))
        availCpus.addAndGet(-result.cpus)
        availMemory.addAndGet(-result.memory)
        return result
    }

    /**
     * Restore the resources associated to the specified reservation
     */
    void release(Reservation reservation) {
        if( isUnlimited() || reservation == null )
            return
        availCpus.addAndGet(reservation.cpus)
        availMemory.addAndGet(reservation.memory)
    }

    /**
     * The amount of resources reserved by a submitted task
     */
    @PackageScope
    static class Reservation {
        final int cpus
        final long memory
        Reservation(int cpus, long memory) {
            this.cpus = cpus
            this.memory = memory
        }
    }

    String toString() {
        isUnlimited()
            ? "ResourceAccount[unlimited]"
            : "ResourceAccount[cpus=${availCpus.get()}/$maxCpus, memory=${new MemoryUnit(availMemory.get())}/${new MemoryUnit(maxMemory)}]"
    }
}
