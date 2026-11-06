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

package nextflow.dataflow.engine

import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.TimeUnit
import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.locks.Condition
import java.util.concurrent.locks.ReentrantLock

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import groovyx.gpars.dataflow.DataflowReadChannel
import groovyx.gpars.dataflow.DataflowWriteChannel
import groovyx.gpars.dataflow.DataflowVariable
import groovyx.gpars.dataflow.operator.DataflowProcessor
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.dag.DAG
import nextflow.extension.CH

/**
 * Bounds the backlog of a run (disabled with {@code NXF_DATAFLOW_HIGH_WATER=0}):
 * - a process stops taking new inputs while its own backlog is at or above the high-water mark;
 * - a source channel blocks its producer thread while any process downstream of it is.
 *
 * Operators never block and operator queues are never counted, so every pause drains
 * without further source input. Only threads flagged with {@link #runAsSource} ever block.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class Backpressure {

    static final private ThreadLocal<Boolean> SOURCE_THREAD = ThreadLocal.withInitial(() -> Boolean.FALSE)

    static final String HIGH_WATER_ENV = 'NXF_DATAFLOW_HIGH_WATER'

    /**
     * The high-water mark, by default twice the dataflow pool size (the most input sets a process
     * prepares at once) and at least 64, so that task preparation alone never pauses a process.
     * Zero disables backpressure
     */
    static int highWater(int poolSize) {
        final value = SysEnv.get(HIGH_WATER_ENV)?.trim()
        if( !value )
            return Math.max(64, 2 * poolSize)
        try {
            final result = Integer.parseInt(value)
            if( result >= 0 )
                return result
        }
        catch( NumberFormatException ignored ) {
        }
        throw new IllegalArgumentException("Invalid $HIGH_WATER_ENV value: '$value' -- it must be a non-negative integer (0 disables backpressure)")
    }

    /**
     * @return The backpressure of the current session, or {@code null} when disabled
     */
    static Backpressure current() {
        final session = Global.session
        return session instanceof Session ? session.getBackpressure() : null
    }

    /**
     * Runs an action with the current thread flagged as a source producer thread
     */
    static void runAsSource(Runnable action) {
        final prev = SOURCE_THREAD.get()
        SOURCE_THREAD.set(Boolean.TRUE)
        try {
            action.run()
        }
        finally {
            SOURCE_THREAD.set(prev)
        }
    }

    /**
     * Blocks a source producer thread while a process fed by the given source channel is over
     * the high-water mark. No-op when disabled or on a thread not flagged with {@link #runAsSource}
     */
    static void awaitSource(Object source) {
        if( SOURCE_THREAD.get() )
            current()?.awaitOpen(source)
    }

    static class Gate {
        final String name
        volatile List<ProcessBacklog> feeds = List.<ProcessBacklog>of()

        // the threads blocked on the gate, and the times a thread blocked on it (guarded by the backpressure lock)
        private volatile int waiters
        private volatile long blockCount

        Gate(String name) { this.name = name }

        boolean isWaiting() { waiters > 0 }

        @PackageScope
        int getWaiters() { waiters }

        @PackageScope
        long getBlockCount() { blockCount }

        @PackageScope
        void enter() {
            this.@waiters = this.@waiters + 1
            this.@blockCount = this.@blockCount + 1
        }

        @PackageScope
        void leave() {
            this.@waiters = this.@waiters - 1
        }

        boolean isOpen(int mark) {
            for( ProcessBacklog it : feeds ) {
                if( it.isOver(mark) )
                    return false
            }
            return true
        }

        String toString() { "$name -> ${feeds*.name}" }
    }

    final int highWater

    final int lowWater

    private final Session session

    private final Map<Object,Gate> gates = Collections.synchronizedMap(new IdentityHashMap<Object,Gate>())

    private final List<ProcessBacklog> processes = new ArrayList<>()

    private final ReentrantLock lock = new ReentrantLock()

    private final Condition changed = lock.newCondition()

    private volatile int waiters

    private volatile boolean released

    private ExecutorService pump

    Backpressure(Session session, int highWater) {
        this.session = session
        this.highWater = highWater
        this.lowWater = highWater.intdiv(2)
        log.debug "Dataflow backpressure enabled -- high-water mark: $highWater; low-water mark: $lowWater"
    }

    boolean isReleased() { released }

    /**
     * Mark a source channel, i.e. one whose producer may be paused. Value channels and topics
     * are never gated, and re-marking a channel keeps its gate
     */
    void markSource(String name, Object source) {
        if( source instanceof Collection ) {
            for( Object it : source ) markSource(name, it)
        }
        else if( source != null && source !instanceof DataflowVariable && !CH.isTopic(source) )
            gates.putIfAbsent(source, new Gate(name))
    }

    @PackageScope
    Gate gate(Object source) { gates.get(source) }

    @PackageScope
    List<Gate> gates() {
        synchronized (gates) {
            return new ArrayList<Gate>(gates.values())
        }
    }

    @PackageScope
    void addProcess(ProcessBacklog backlog) {
        synchronized (processes) {
            processes.add(backlog)
        }
    }

    @PackageScope
    void removeProcess(ProcessBacklog backlog) {
        synchronized (processes) {
            processes.remove(backlog)
        }
    }

    ProcessBacklog register(String name, List inputs, DataflowProcessor operator) {
        // the last input is the control channel
        final queues = new ArrayList<DataflowReadChannel>(inputs.subList(0, inputs.size()-1) as List<DataflowReadChannel>)
        final result = new ProcessBacklog(this, name, queues, operator)
        operator.setGate(result)
        addProcess(result)
        return result
    }

    /**
     * Compute the processes downstream of each source
     */
    void connect(DAG dag) {
        final next = new HashMap<DAG.Vertex,List<DAG.Vertex>>()
        final origin = new IdentityHashMap<Object,DAG.Vertex>()
        for( DAG.Edge edge : dag.edges ) {
            if( edge.from == null )
                continue
            origin.put(edge.channel, edge.from)
            if( edge.to != null )
                next.computeIfAbsent(edge.from, (k) -> new ArrayList<DAG.Vertex>()).add(edge.to)
        }
        synchronized (gates) {
            for( Map.Entry<Object,Gate> entry : gates.entrySet() ) {
                final start = origin.get(entry.key)
                if( start != null )
                    entry.value.feeds = feeds(start, next)
                log.debug "Backpressure source ${entry.value}"
            }
        }
    }

    static private List<ProcessBacklog> feeds(DAG.Vertex start, Map<DAG.Vertex,List<DAG.Vertex>> next) {
        final result = new ArrayList<ProcessBacklog>()
        final visited = new HashSet<DAG.Vertex>()
        final queue = new ArrayDeque<DAG.Vertex>()
        queue.add(start)
        while( !queue.isEmpty() ) {
            final vertex = queue.poll()
            if( !visited.add(vertex) )
                continue
            final backlog = vertex.process?.getBacklog()
            if( backlog != null )
                result.add(backlog)
            final targets = next.get(vertex)
            if( targets )
                queue.addAll(targets)
        }
        return result
    }

    boolean isOpen(Object source) {
        final gate = gates.get(source)
        return gate == null || released || gate.isOpen(highWater)
    }

    /**
     * Block the current thread when a process fed by the given source is at the high-water mark,
     * until every one is below the low-water mark. No-op on a thread not flagged as a source thread
     */
    void awaitOpen(Object source) {
        if( !SOURCE_THREAD.get() )
            return
        final gate = gates.get(source)
        if( gate == null || released || gate.isOpen(highWater) )
            return
        lock.lock()
        try {
            waiters++
            gate.enter()
            while( !released && !gate.isOpen(lowWater+1) )
                changed.await(1, TimeUnit.SECONDS)
        }
        catch( InterruptedException e ) {
            Thread.currentThread().interrupt()
        }
        finally {
            gate.leave()
            waiters--
            lock.unlock()
        }
    }

    /**
     * Emit the given items into a source channel. At most a high-water mark of items is bound
     * on the calling thread while the gate is open, then the rest is handed to a source thread,
     * so the caller never blocks. Items bound on the calling thread are not paced by operators
     * downstream, hence the limit
     */
    void emit(DataflowWriteChannel ch, Iterable items, boolean close) {
        // like a Groovy for loop, a null iterable emits nothing
        final itr = items != null ? items.iterator() : Collections.emptyIterator()
        final gate = gates.get(ch)
        final limit = gate?.feeds ? highWater : Integer.MAX_VALUE
        int count = 0
        while( itr.hasNext() && count++ < limit && isOpen(ch) )
            ch.bind(itr.next())
        if( !itr.hasNext() ) {
            if( close ) ch.bind(CH.stop())
            return
        }
        log.trace "Source ${gate} -- emitting remaining items on a source thread"
        pump().execute(() -> runAsSource(() -> emit0(ch, itr, close)))
    }

    private void emit0(DataflowWriteChannel ch, Iterator itr, boolean close) {
        try {
            while( itr.hasNext() ) {
                awaitOpen(ch)
                ch.bind(itr.next())
            }
            if( close ) ch.bind(CH.stop())
        }
        catch( Throwable e ) {
            log.error("Unexpected error while emitting source items", e)
            session.abort(e)
        }
    }

    private synchronized ExecutorService pump() {
        if( pump == null ) {
            final count = new AtomicInteger()
            pump = Executors.newCachedThreadPool((Runnable r) -> {
                final thread = new Thread(r, "SourcePump-${count.incrementAndGet()}")
                thread.setDaemon(true)
                return thread
            })
        }
        return pump
    }

    /**
     * Wake the sources waiting on a gate, to be called after a backlog decrement
     */
    void wake() {
        if( waiters > 0 )
            signalWaiters()
    }

    private void signalWaiters() {
        lock.lock()
        try {
            changed.signalAll()
        }
        finally {
            lock.unlock()
        }
    }

    /**
     * Open every gate for the rest of the run, e.g. on abort or cancel
     */
    void release() {
        if( released )
            return
        released = true
        for( ProcessBacklog it : processes() )
            it.operator?.resume()
        signalWaiters()
    }

    void shutdown() {
        release()
        synchronized (this) {
            pump?.shutdown()
        }
    }

    @PackageScope
    List<ProcessBacklog> processes() {
        synchronized (processes) {
            return new ArrayList<ProcessBacklog>(processes)
        }
    }
}
