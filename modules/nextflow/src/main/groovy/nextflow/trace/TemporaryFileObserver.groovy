/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.trace

import java.nio.file.Path
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.locks.Lock
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.dag.DAG
import nextflow.file.FileHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.params.FileOutParam
/**
 * Delete task directories once they are no longer needed.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TemporaryFileObserver implements TraceObserver {

    private DAG dag

    private Map<String,ProcessState> processes = new HashMap<>()

    private Map<TaskRun,TaskState> tasks = new HashMap<>()

    private Map<Path,TaskRun> taskLookup = new HashMap<>()

    private Lock sync = new ReentrantLock()

    @Override
    void onFlowCreate(Session session) {
        this.dag = session.dag
    }

    /**
     * When the workflow begins, determine the consumers of each process
     * in the DAG.
     */
    @Override
    void onFlowBegin() {

        for( def processNode : dag.vertices ) {
            // skip nodes that are not processes
            if( !processNode.process )
                continue

            // find all downstream processes in the abstract dag
            def processName = processNode.process.name
            def consumers = [] as Set
            def queue = [ processNode ]

            while( !queue.isEmpty() ) {
                // remove a node from the search queue
                final sourceNode = queue.remove(0)

                // search each outgoing edge from the source node
                for( def edge : dag.edges ) {
                    if( edge.from != sourceNode )
                        continue

                    def node = edge.to

                    // skip if process is terminal
                    if( !node )
                        continue

                    // add process nodes to the list of consumers
                    if( node.process != null )
                        consumers.add(node.process.name)
                    // add operator nodes to the queue to keep searching
                    else
                        queue.add(node)
                }
            }

            log.trace "Process `${processName}` is consumed by the following processes: ${consumers.collect({ "`${it}`" }).join(', ')}"

            processes[processName] = new ProcessState(consumers ?: [processName] as Set)
        }
    }

    /**
     * When a task is created, add it to the state map and add it as a consumer
     * of any task whose output it takes as input.
     *
     * @param handler
     * @param trace
     */
    @Override
    void onProcessPending(TaskHandler handler, TraceRecord trace) {
        // query task input files
        final task = handler.task
        final inputs = task.getInputFilesMap().values()

        sync.lock()
        try {
            // add task to the task state map
            tasks[task] = new TaskState()

            // add task as consumer to each task whoes output it takes as input
            for( Path path : inputs )
                if( path in taskLookup )
                    tasks[taskLookup[path]].consumers.add(task)
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * When a task is completed, track any temporary output files
     * for automatic cleanup.
     *
     * @param handler
     * @param trace
     */
    @Override
    void onProcessComplete(TaskHandler handler, TraceRecord trace) {
        // query task output files
        final task = handler.task
        final outputs = task
            .getOutputsByType(FileOutParam)
            .values()
            .flatten() as List<Path>

        sync.lock()
        try {
            // mark task as completed
            tasks[task].completed = true

            // scan tasks for cleanup
            cleanup0()

            // add new output files to task lookup
            for( Path path : outputs )
                taskLookup[path] = task
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * When a process is closed (all tasks of the process have been created),
     * mark the process as closed and scan tasks for cleanup.
     *
     * @param process
     */
    @Override
    void onProcessClose(TaskProcessor process) {
        sync.lock()
        try {
            processes[process.name].closed = true
            cleanup0()
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Delete any task directories that have no more barriers.
     */
    private void cleanup0() {
        for( TaskRun task : tasks.keySet() ) {
            final taskState = tasks[task]
            if( taskState.completed && !taskState.deleted && canDelete(task) ) {
                log.trace "Deleting task directory: ${task.workDir.toUriString()}"
                FileHelper.deletePath(task.workDir)
                taskState.deleted = true
            }
        }
    }

    /**
     * Determine whether a task directory can be deleted.
     *
     * A task directory can be deleted if all of its process consumers
     * are closed and all of its task consumers are completed.
     */
    private boolean canDelete(TaskRun task) {
        final taskState = tasks[task]
        final processConsumers = processes[task.processor.name].consumers
        final taskConsumers = taskState.consumers
        processConsumers.every { p -> processes[p].closed } && taskConsumers.every { t -> tasks[t].completed }
    }

    static private class ProcessState {
        Set<String> consumers
        boolean closed = false

        ProcessState(Set<String> consumers) {
            this.consumers = consumers
        }
    }

    static private class TaskState {
        Set<TaskRun> consumers = [] as Set
        boolean completed = false
        boolean deleted = false
    }

}
