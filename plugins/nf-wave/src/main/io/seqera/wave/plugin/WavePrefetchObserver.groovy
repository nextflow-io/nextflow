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
import java.util.concurrent.LinkedBlockingQueue
import java.util.concurrent.ThreadPoolExecutor
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.dag.DAG
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.trace.TraceObserverV2
import nextflow.trace.event.TaskEvent
import nextflow.util.Duration

/**
 * Request ahead the containers of the processes downstream of a running process, so that their
 * Wave provisioning (e.g. the SOCI index of the image, or a container build) runs while the
 * upstream tasks run, rather than when their own first task is submitted.
 *
 * When the first task of a process is pending, the processes up to {@code depth} levels
 * downstream in the DAG are resolved with a task preview, the same way as {@code nextflow inspect}.
 * The requests go through a small pool of threads, and each process is requested at most once.
 * As with {@code nextflow inspect -concretize}, a thread awaits the container to be ready before
 * taking the next request, so that at most {@code maxConcurrent} containers are provisioned
 * ahead at the same time, the nearest ones first. A queued process whose own task starts in the
 * meantime is skipped.
 * The Wave client caches the response, so the actual task reuses it. A process whose container
 * depends on the task inputs cannot be resolved ahead, and it is requested by its task as usual.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class WavePrefetchObserver implements TraceObserverV2 {

    /** default max number of containers provisioned ahead at the same time */
    static final int DEFAULT_MAX_CONCURRENT = 4

    private DAG dag

    private int depth

    /** max number of containers provisioned ahead at the same time */
    private int maxConcurrent

    private ThreadPoolExecutor executor

    /** the processes whose downstream was already looked up */
    private Set<String> triggered = ConcurrentHashMap.newKeySet()

    /** the processes whose container was already requested, ahead or by their own task */
    private Set<String> requested = ConcurrentHashMap.newKeySet()

    /** the processes having a task pending, therefore requesting the container on their own */
    private Set<String> started = ConcurrentHashMap.newKeySet()

    /** the process whose first task caused a process to be requested ahead, for debugging */
    private Map<String,String> triggers = new ConcurrentHashMap<>()

    private AtomicInteger countRequested = new AtomicInteger()
    private AtomicInteger countReady = new AtomicInteger()
    private AtomicInteger countFailed = new AtomicInteger()
    private AtomicInteger countSkipped = new AtomicInteger()

    WavePrefetchObserver(DAG dag, int depth) {
        this(dag, depth, DEFAULT_MAX_CONCURRENT)
    }

    WavePrefetchObserver(DAG dag, int depth, int maxConcurrent) {
        this.dag = dag
        this.depth = depth
        this.maxConcurrent = maxConcurrent
        this.executor = new ThreadPoolExecutor(maxConcurrent, maxConcurrent, 0L, TimeUnit.MILLISECONDS, new LinkedBlockingQueue<Runnable>(), (Runnable r) -> {
            final thread = new Thread(r, 'wave-prefetch')
            thread.setDaemon(true)
            return thread
        })
    }

    @Override
    void onTaskPending(TaskEvent event) {
        final process = event.handler.task.processor
        // its own task requests the container of this process
        started.add(process.name)
        requested.add(process.name)
        if( !triggered.add(process.name) )
            return
        final downstream = downstream(dag, process, depth)
        final submitted = new ArrayList<String>()
        for( TaskProcessor next : downstream ) {
            if( requested.add(next.name) ) {
                triggers.put(next.name, process.name)
                submitted.add(next.name)
                submit(next)
            }
        }
        if( log.isDebugEnabled() ) {
            log.debug "Wave prefetch triggered by process '${process.name}' - downstream (depth ${depth}): ${downstream*.name}; new to request: ${submitted}; ${slots()}"
        }
    }

    protected void submit(TaskProcessor process) {
        final queuedAt = System.currentTimeMillis()
        executor.submit(() -> prefetch(process, queuedAt))
    }

    @Override
    void onFlowComplete() {
        log.debug "Wave prefetch completed - triggered by: ${triggered.size()} processes; requested ahead: ${countRequested}; ready: ${countReady}; failed: ${countFailed}; skipped: ${countSkipped}; ${slots()}"
        // interrupt the threads still awaiting a container
        executor.shutdownNow()
    }

    protected void prefetch(TaskProcessor process) {
        prefetch(process, System.currentTimeMillis())
    }

    protected void prefetch(TaskProcessor process, long queuedAt) {
        final trigger = triggers.get(process.name)
        final startedAt = System.currentTimeMillis()
        if( started.contains(process.name) ) {
            countSkipped.incrementAndGet()
            log.debug "Skip requesting ahead the Wave container of process '${process.name}' - its task already started (trigger: '${trigger}'; queued: ${elapsed(queuedAt, startedAt)})"
            return
        }
        try {
            final task = process.createTaskPreview()
            final container = task.getContainer()
            final requestedAt = System.currentTimeMillis()
            countRequested.incrementAndGet()
            log.debug "Wave container requested ahead for process '${process.name}': ${container} (trigger: '${trigger}'; queued: ${elapsed(queuedAt, startedAt)}; request: ${elapsed(startedAt, requestedAt)}; ${describe(task, startedAt)}; ${slots()})"
            // keep the thread busy until the container is ready, so the next request waits for it
            final checks = awaitReady(task)
            countReady.incrementAndGet()
            log.debug "Wave container ready ahead for process '${process.name}': ${container} (await: ${elapsed(requestedAt, System.currentTimeMillis())}; status checks: ${checks})"
        }
        catch( InterruptedException e ) {
            Thread.currentThread().interrupt()
        }
        catch( Throwable e ) {
            // e.g. a container depending on the task inputs, the task requests it as usual
            countFailed.incrementAndGet()
            log.debug "Unable to request ahead the Wave container of process '${process.name}' (trigger: '${trigger}'; elapsed: ${elapsed(startedAt, System.currentTimeMillis())}) - cause: ${e.message}"
        }
    }

    /**
     * Await the task container to be ready, polling it the same way as {@code nextflow inspect -concretize}.
     * A failed or timed out provisioning is thrown by the Wave client.
     *
     * @return The number of times the container status was checked
     */
    protected int awaitReady(TaskRun task) {
        int c=0
        while( !task.isContainerReady() ) {
            sleep(c++ < 10 ? 1_000 : 5_000)
        }
        return c+1
    }

    /**
     * The Wave metadata of the task container, to match it with the Wave service logs
     */
    @PackageScope
    static String describe(TaskRun task, long startedAt) {
        final meta = task.containerMeta()
        if( !meta )
            return 'meta: none'
        // the Wave client caches the response by the request assets, so a container already requested is not sent again
        final reused = meta.requestTime && meta.requestTime.toInstant().toEpochMilli() < startedAt
        return "response: ${reused ? 'reused' : 'new'}; requestId: ${meta.requestId}; buildId: ${meta.buildId ?: meta.mirrorId}; cached: ${meta.cached}; scanId: ${meta.scanId}; source: ${meta.sourceImage}"
    }

    private String slots() {
        return "slots running: ${executor.getActiveCount()}/${maxConcurrent}; queued: ${executor.getQueue().size()}"
    }

    static private String elapsed(long from, long to) {
        return Duration.of(to-from).toString()
    }

    protected void sleep(long millis) {
        Thread.sleep(millis)
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
