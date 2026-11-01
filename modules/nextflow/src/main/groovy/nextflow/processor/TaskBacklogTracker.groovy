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

import java.util.concurrent.ConcurrentHashMap
import java.util.function.Supplier

import groovy.transform.CompileStatic
import nextflow.dataflow.engine.ProcessBacklog

/**
 * Counts the tasks of a process in its backlog, from the input set taken by the
 * process to the submission of its task by the executor:
 *
 * <pre>
 * accepted --> preparing --> scheduled --> pending --> dequeued
 *                  \
 *                   +--> prepared (discarded: `when` guard, cached, stored, failed, array child)
 * </pre>
 *
 * It is the only component that changes the task counts of a {@link ProcessBacklog}.
 * Each task is tracked by identity and every transition is idempotent, so a task is
 * counted at most once whatever the path, and a task that was never counted (e.g. a
 * task submitted without backpressure) is ignored.
 *
 * A retried task and a task array are scheduled without being prepared, so they are
 * counted as pending only.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class TaskBacklogTracker {

    private final Supplier<ProcessBacklog> backlog

    private final Set<TaskRun> preparing = ConcurrentHashMap.newKeySet()

    private final Set<TaskRun> pending = ConcurrentHashMap.newKeySet()

    TaskBacklogTracker(Supplier<ProcessBacklog> backlog) {
        this.backlog = backlog
    }

    /**
     * Count a task as preparing, from the moment the process takes its input set
     */
    void accepted(TaskRun task) {
        final backlog = backlog.get()
        if( backlog == null || !preparing.add(task) )
            return
        backlog.taskStarted()
    }

    /**
     * Count a task as pending, before it is handed to the executor. A preparing task
     * moves from preparing to pending at once
     */
    void scheduled(TaskRun task) {
        final backlog = backlog.get()
        if( backlog == null || !pending.add(task) )
            return
        if( preparing.remove(task) )
            backlog.taskScheduled()
        else
            backlog.taskPending()
    }

    /**
     * End the preparation of a task: release it unless it was scheduled
     */
    void prepared(TaskRun task) {
        if( preparing.remove(task) )
            backlog.get()?.taskPrepared()
    }

    /**
     * Release a pending task, once the executor takes it off its pending queue
     * (either submitted or failed), or when the hand-over to the executor fails
     */
    void dequeued(TaskRun task) {
        if( pending.remove(task) )
            backlog.get()?.taskSubmitted()
    }

    int getPreparing() { preparing.size() }

    int getPending() { pending.size() }
}
