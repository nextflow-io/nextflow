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
import nextflow.dag.ConcreteDAG
import nextflow.file.FileHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.params.FileOutParam
/**
 * Track temporary output files and delete them once they
 * are no longer needed.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class TemporaryFileObserver implements TraceObserver {

    private DAG dag

    private Map<Path,Status> pathStatuses = new HashMap<>()

    private Set<TaskRun> completedTasks = new HashSet<>()

    private Map<String,Set<String>> processConsumers = new HashMap<>()

    private Set<String> closedProcesses = new HashSet<>()

    private Lock sync = new ReentrantLock()

    @Override
    void onFlowCreate(Session session) {
        this.dag = session.dag
    }

    /**
     * When a task is created, add it to the set of consumers for any temporary
     * file that it takes as input.
     *
     * @param handler
     * @param trace
     */
    @Override
    void onProcessPending(TaskHandler handler, TraceRecord trace) {
        sync.lock()
        try {
            final task = handler.task
            for( def entry : task.getInputFilesMap() ) {
                final storePath = entry.value
                if( storePath in pathStatuses )
                    pathStatuses[storePath].consumers.add(task)
            }
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
        // query all temporary files produced by task
        final task = handler.task
        final tempOutputs = task
            .getOutputsByType(FileOutParam)
            .findAll { param, paths -> param.temporary }
            .values()
            .flatten() as Set<Path>

        sync.lock()
        try {
            // mark task as completed
            completedTasks.add(task)

            // scan temporary files for cleanup
            cleanup0()

            // add new temporary outputs to status map
            for( Path path : tempOutputs )
                pathStatuses[path] = new Status(task)
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Get the consumers of a process.
     *
     * @param processName
     */
    private Set<String> getProcessConsumers(String processName) {
        if( processName !in processConsumers )
            processConsumers[processName] = getProcessConsumers0(processName)

        return processConsumers[processName]
    }

    private Set<String> getProcessConsumers0(String processName) {

        // find the task's process node in the abstract dag
        final processNode = dag.vertices
            .find { node -> node.process?.name == processName }

        // find all downstream processes in the abstract dag
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

                // add process nodes to the list of consumers
                if( node.process != null )
                    consumers.add(node.process.name)
                // add operator nodes to the queue to keep searching
                else
                    queue.add(node)
            }
        }

        log.trace "Process ${processName} has the following consumers: ${consumers.join(', ')}"

        return consumers
    }

    /**
     * When a process is closed (all tasks of the process have been created),
     * mark the process as closed and scan for automatic cleanup.
     *
     * @param process
     */
    @Override
    void onProcessClose(TaskProcessor process) {
        sync.lock()
        try {
            closedProcesses.add(process.name)
            cleanup0()
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Delete any temporary file that has no more barriers.
     */
    private void cleanup0() {
        for( Path path : pathStatuses.keySet() ) {
            final status = pathStatuses[path]
            if( !status.deleted && canDelete(path) ) {
                log.trace "Deleting temporary file: ${path}"
                FileHelper.deletePath(path)
                status.deleted = true
            }
        }
    }

    /**
     * Determine whether a path can be deleted.
     *
     * A path can be deleted if all of its process consumers
     * are closed and all of its task consumers are completed.
     */
    private boolean canDelete(Path path) {
        final status = pathStatuses[path]
        final processConsumers = getProcessConsumers(status.task.processor.name)
        final taskConsumers = status.consumers
        processConsumers.every { p -> p in closedProcesses } && taskConsumers.every { t -> t in completedTasks }
    }

    static private class Status {
        TaskRun task
        Set<TaskRun> consumers = [] as Set
        boolean deleted = false

        Status(TaskRun task) {
            this.task = task
        }
    }

}
