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

import groovy.transform.MapConstructor
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.dag.DAG
import nextflow.file.FileHelper
import nextflow.processor.TaskHandler
import nextflow.processor.TaskProcessor
/**
 * Watch for temporary output files and "clean" them once they
 * are no longer needed.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
class TemporaryFileObserver implements TraceObserver {

    private DAG dag

    private Map<String,Status> statusMap

    @Override
    void onFlowCreate(Session session) {
        this.dag = session.dag
        this.statusMap = new ConcurrentHashMap<>()
    }

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

            statusMap[processName] = new Status(
                paths: [] as Set,
                processBarriers: consumers ?: [processName] as Set
            )
        }
    }

    /**
     * When a task is started, track the task directory for automatic cleanup.
     *
     * @param handler
     * @param trace
     */
    @Override
    synchronized void onProcessStart(TaskHandler handler, TraceRecord trace) {
        // add task directory to status map
        final task = handler.task
        final processName = task.processor.name

        log.trace "Task completed from process `${processName}`, tracking task directory for automatic cleanup"

        statusMap[processName].paths.add(task.workDir)
    }

    /**
     * When all tasks of a process are completed, update the status of any processes
     * that are waiting on it in order to be cleaned up.
     *
     * If, after this update is removed, a process has no more barriers,
     * then delete all temporary files for that process.
     *
     * @param process
     */
    @Override
    synchronized void onProcessTerminate(TaskProcessor process) {
        log.trace "Process `${process.name}` is complete, updating barriers for upstream processes"

        for( def entry : statusMap ) {
            // remove barrier from each upstream process
            final producer = entry.key
            final status = entry.value
            final barriers = status.processBarriers

            barriers.remove(process.name)

            // if a process has no more barriers, trigger the cleanup
            if( barriers.isEmpty() ) {
                log.trace "All barriers for process `${producer}` are complete, time to clean up temporary files"

                deleteTemporaryFiles(status.paths)
                statusMap.remove(producer)
            }
        }
    }

    void deleteTemporaryFiles(Collection<Path> paths) {
        for( Path path : paths ) {
            log.trace "Cleaning temporary file: ${path}"

            FileHelper.deletePath(path)
        }
    }

    @MapConstructor
    static protected class Status {
        Set<Path> paths
        Set<String> processBarriers
    }

}
