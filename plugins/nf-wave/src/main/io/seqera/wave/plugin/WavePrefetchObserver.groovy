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

package io.seqera.wave.plugin

import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.dag.DAG
import nextflow.processor.TaskProcessor
import nextflow.trace.TraceObserverV2
import nextflow.trace.event.TaskEvent

/**
 * Request ahead the containers of the processes downstream of a running process, so that their
 * Wave provisioning (e.g. the SOCI index of the image, or a container build) runs while the
 * upstream tasks run, rather than when their own first task is submitted.
 *
 * When the first task of a process is pending, the processes up to {@code depth} levels
 * downstream in the DAG are resolved with a task preview, the same way as {@code nextflow inspect}.
 * The requests go through a small pool of threads, and each process is requested at most once.
 * The Wave client caches the response, so the actual task reuses it. A process whose container
 * depends on the task inputs cannot be resolved ahead, and it is requested by its task as usual.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WavePrefetchObserver implements TraceObserverV2 {

    /** max number of container requests sent ahead at the same time */
    static final int MAX_CONCURRENT = 4

    private DAG dag

    private int depth

    private ExecutorService executor

    /** the processes whose downstream was already looked up */
    private Set<String> triggered = ConcurrentHashMap.newKeySet()

    /** the processes whose container was already requested, ahead or by their own task */
    private Set<String> requested = ConcurrentHashMap.newKeySet()

    WavePrefetchObserver(DAG dag, int depth) {
        this.dag = dag
        this.depth = depth
        this.executor = Executors.newFixedThreadPool(MAX_CONCURRENT, (Runnable r) -> {
            final thread = new Thread(r, 'wave-prefetch')
            thread.setDaemon(true)
            return thread
        })
    }

    @Override
    void onTaskPending(TaskEvent event) {
        final process = event.handler.task.processor
        // its own task requests the container of this process
        requested.add(process.name)
        if( !triggered.add(process.name) )
            return
        for( TaskProcessor next : downstream(dag, process, depth) ) {
            if( requested.add(next.name) )
                submit(next)
        }
    }

    protected void submit(TaskProcessor process) {
        executor.submit(() -> prefetch(process))
    }

    @Override
    void onFlowComplete() {
        executor.shutdown()
    }

    protected void prefetch(TaskProcessor process) {
        try {
            final container = process.createTaskPreview().getContainer()
            log.debug "Wave container requested ahead for process '${process.name}': ${container}"
        }
        catch( Throwable e ) {
            // e.g. a container depending on the task inputs, the task requests it as usual
            log.debug "Unable to request ahead the Wave container of process '${process.name}' - cause: ${e.message}"
        }
    }

    /**
     * The processes up to {@code depth} levels downstream of the given process. Operators between
     * two processes do not count as a level.
     *
     * @param dag The workflow DAG
     * @param process The process from which the lookup starts
     * @param depth The number of process levels to look up
     * @return The downstream processes, nearest first
     */
    @PackageScope
    static List<TaskProcessor> downstream(DAG dag, TaskProcessor process, int depth) {
        final result = new LinkedHashSet<TaskProcessor>()
        final successors = dag.edges.groupBy { DAG.Edge e -> e.from }
        List<DAG.Vertex> frontier = dag.vertices.findAll { DAG.Vertex v -> v.process.is(process) }
        for( int level=0; level<depth && frontier; level++ ) {
            final next = new ArrayList<DAG.Vertex>()
            for( DAG.Vertex v : frontier ) {
                for( DAG.Vertex p : nextProcesses(successors, v) ) {
                    if( !p.process.is(process) && result.add(p.process) )
                        next.add(p)
                }
            }
            frontier = next
        }
        return new ArrayList<TaskProcessor>(result)
    }

    /**
     * The process vertices fed by the given vertex, through any number of non-process vertices
     */
    static private List<DAG.Vertex> nextProcesses(Map<DAG.Vertex,List<DAG.Edge>> successors, DAG.Vertex from) {
        final result = new ArrayList<DAG.Vertex>()
        final visited = new HashSet<DAG.Vertex>()
        final queue = new LinkedList<DAG.Vertex>([from])
        while( queue ) {
            final v = queue.removeFirst()
            for( DAG.Edge e : successors.getOrDefault(v, List.<DAG.Edge>of()) ) {
                if( !e.to || !visited.add(e.to) )
                    continue
                // a process vertex holds its processor, see NodeMarker
                if( e.to.process != null )
                    result.add(e.to)
                else
                    queue.add(e.to)
            }
        }
        return result
    }
}
