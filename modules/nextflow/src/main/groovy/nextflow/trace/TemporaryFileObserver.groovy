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

import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.dag.DAG
import nextflow.dag.ConcreteDAG
import nextflow.file.FileHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
import nextflow.script.params.FileOutParam
/**
 * Track temporary output files and delete them once they
 * are no longer needed.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
class TemporaryFileObserver implements TraceObserver {

    private DAG dag

    private Map<String,Status> statusMap = new ConcurrentHashMap<>()

    private Lock sync = new ReentrantLock()

    @Override
    void onFlowCreate(Session session) {
        this.dag = session.dag
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
        // find all temporary output files
        final task = handler.task
        final tempOutputs = task
            .getOutputsByType(FileOutParam)
            .findAll { param, paths -> param.temporary }
            .values()
            .flatten()

        if( tempOutputs.isEmpty() )
            return

        // update status tracker for the task's process
        final processName = task.processor.name

        sync.lock()
        try {
            if( processName !in statusMap ) {
                log.trace "Process ${processName} has temporary output files, tracking for automatic cleanup"
                statusMap[processName] = new Status(findAllConsumers(processName))
            }
            statusMap[processName].paths.addAll(tempOutputs)
        }
        finally {
            sync.unlock()
        }
    }

    /**
     * Find all processes which are consumers of a given process.
     *
     * @param processName
     */
    private Set<String> findAllConsumers(String processName) {

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
     * When a process is completed, update the status of any processes
     * that are waiting on it in order to cleanup temporary outputs.
     *
     * If, after this update is removed, a process has no more barriers,
     * then clean all temporary files for that process.
     *
     * @param process
     */
    @Override
    void onProcessTerminate(TaskProcessor process) {
        sync.lock()
        try {
            for( def entry : statusMap ) {
                // remove barrier from each upstream process
                final producer = entry.key
                final status = entry.value
                final consumers = status.processBarriers

                consumers.remove(process.name)

                // if a process has no more barriers, trigger the cleanup
                if( consumers.isEmpty() ) {
                    log.trace "All consumers of process ${producer} are complete, deleting temporary files"

                    deleteTemporaryFiles(status.paths)
                    statusMap.remove(producer)
                }
            }
        }
        finally {
            sync.unlock()
        }
    }

    private void deleteTemporaryFiles(Collection<Path> paths) {
        for( Path path : paths ) {
            log.trace "Deleting temporary file: ${path}"
            FileHelper.deletePath(path)
        }
    }

    static private class Status {
        Set<Path> paths
        Set<String> processBarriers

        Status(Set<String> processes) {
            this.paths = [] as Set
            this.processBarriers = processes
        }
    }

}
