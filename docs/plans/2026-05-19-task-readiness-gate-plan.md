# TaskReadinessGate Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Add a `TaskReadinessGate` plugin extension point to Nextflow core that defers task submission until external preconditions are met (e.g. Glacier restore), running gate work on a managed virtual-thread executor.

**Architecture:** Plugins implement a `void prepare(TaskHandler)` method that may block freely. `TaskPollingMonitor` submits each gate's `prepare` to a managed `ExecutorService` when the task is scheduled, polls the resulting `Future` from inside `canSubmit`, and routes any thrown exception through the existing `handleException` failure path so `errorStrategy` applies. An `executor.gateMaxWait` config option (default `24h`) bounds wait time.

**Tech Stack:** Groovy 4.0.29, Spock framework, PF4J plugin loading (`Plugins.getExtensions`), `java.util.concurrent` (virtual threads via `Threads.useVirtual()`).

**Reference:** ADR `adr/20260516-task-readiness-gate.md`, PR [#7151](https://github.com/nextflow-io/nextflow/pull/7151), implementation design `docs/plans/2026-05-19-task-readiness-gate-implementation.md`.

**Working directory:** `/Users/robsyme/dev/bears/nextflow/master/.worktrees/task-readiness-gate`. All paths below are relative to this directory.

**Conventions:**
- Every commit uses `git commit -s` (DCO sign-off required upstream).
- Apache 2.0 license headers on all new files (copy from any existing source file in the same module).
- Spock specs in `modules/nextflow/src/test/groovy/...` mirror the main source tree.
- `@CompileStatic` on new classes where the parent contract allows it.

---

## Task 1: Add the `TaskReadinessGate` SPI interface

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/TaskReadinessGate.groovy`
- Create: `modules/nextflow/src/test/groovy/nextflow/processor/TaskReadinessGateTest.groovy`

**Step 1: Write the failing test**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * ...license header...
 */
package nextflow.processor

import org.pf4j.ExtensionPoint
import spock.lang.Specification

class TaskReadinessGateTest extends Specification {

    def 'should be a pf4j ExtensionPoint'() {
        expect:
        ExtensionPoint.isAssignableFrom(TaskReadinessGate)
    }

    def 'should declare prepare(TaskHandler) returning void'() {
        when:
        def method = TaskReadinessGate.getMethod('prepare', TaskHandler)
        then:
        method.returnType == void.class
    }
}
```

The canonical extension-point base for Nextflow plugins is `org.pf4j.ExtensionPoint` (every other `extends ExtensionPoint` declaration in the codebase imports it from there — see `nextflow.trace.TraceObserverFactory`, `nextflow.cache.CacheFactory`, etc.).

**Step 2: Run test to verify it fails**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskReadinessGateTest'`
Expected: compilation error — `TaskReadinessGate` does not exist.

**Step 3: Write minimal implementation**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * ...license header...
 */
package nextflow.processor

import groovy.transform.CompileStatic
import org.pf4j.ExtensionPoint

/**
 * Plugin extension point that defers task submission until an external precondition is met.
 *
 * <p>Implementations are invoked once per task by the scheduler, on a managed background
 * thread (virtual thread when available). The task is admitted for submission when this
 * method returns; throw to mark the task as permanently failed and route the cause through
 * the task's {@code errorStrategy}.
 *
 * <p>Implementations may block freely — {@code Thread.sleep}, network calls, long polling.
 * The scheduler thread is never blocked by this call.
 *
 * <p>Implementations <b>must</b> honor {@code Thread.interrupt()} so that task eviction,
 * workflow abort, and the {@code executor.gateMaxWait} backstop can unblock {@code prepare}
 * promptly. Use interruptible primitives ({@code Thread.sleep}, blocking I/O on NIO
 * channels, {@code Future.get}) and propagate {@code InterruptedException}.
 *
 * <p>When multiple gates are registered, a task is admitted only when every gate's
 * {@code prepare} method has returned successfully. Evaluation order across gates is
 * unspecified; all gates start in parallel on the managed executor.
 *
 * <p>Per-process opt-out is available via the {@code hints} directive — gates that wish
 * to support it should check a namespaced hint key (e.g. {@code 'glacier/skip': true}) and
 * return immediately when set. No core mechanism is required.
 */
@CompileStatic
interface TaskReadinessGate extends ExtensionPoint {
    void prepare(TaskHandler handler) throws InterruptedException
}
```

**Step 4: Run test to verify it passes**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskReadinessGateTest'`
Expected: 2 tests pass.

**Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/TaskReadinessGate.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/TaskReadinessGateTest.groovy
git commit -s -m "Add TaskReadinessGate plugin extension point interface"
```

---

## Task 2: Add `executor.gateMaxWait` config option

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/executor/ExecutorConfig.groovy`
- Modify: `modules/nextflow/src/test/groovy/nextflow/executor/ExecutorConfigTest.groovy` (if it exists; otherwise create)

**Step 1: Write the failing test**

Locate the existing test (likely present given the pattern of other config classes). Add:

```groovy
def 'should default gateMaxWait to 24h'() {
    when:
    def config = new ExecutorConfig([:])
    then:
    config.gateMaxWait == Duration.of('24h')
}

def 'should override gateMaxWait from opts'() {
    when:
    def config = new ExecutorConfig(gateMaxWait: '48h')
    then:
    config.gateMaxWait == Duration.of('48h')
}
```

If no test file exists, create one mirroring the style of `TaskPollingMonitorTest`. (Note: there is no "wait indefinitely" sentinel — the standard `opts.foo as Duration ?: default` idiom used throughout `ExecutorConfig` collapses null into the default. If a wait-forever option is later wanted, it would need either a separate boolean or a sentinel Duration value; out of scope for v1.)

**Step 2: Run test to verify it fails**

Run: `./gradlew :nextflow:test --tests 'nextflow.executor.ExecutorConfigTest'`
Expected: FAIL — `gateMaxWait` property does not exist.

**Step 3: Implement**

In `ExecutorConfig.groovy`, alphabetically slot a new field between `exitReadTimeout` and `jobName` (declaration only — no field initialiser):

```groovy
    @ConfigOption
    @Description("""
        Maximum time a `TaskReadinessGate` plugin may take to prepare a task before the task is failed (default: `24h`).
    """)
    final Duration gateMaxWait
```

Then wire the default in the `ExecutorConfig(Map opts)` constructor, alphabetically next to the existing `Duration` defaults at lines 174-187:

```groovy
        gateMaxWait = opts.gateMaxWait as Duration ?: Duration.of('24h')
```

This matches the established pattern used by `dumpInterval`, `exitReadTimeout`, `queueStatInterval`, etc.

**Step 4: Run test to verify it passes**

Run: `./gradlew :nextflow:test --tests 'nextflow.executor.ExecutorConfigTest'`
Expected: 3 new tests pass; existing tests still pass.

**Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/executor/ExecutorConfig.groovy \
        modules/nextflow/src/test/groovy/nextflow/executor/ExecutorConfigTest.groovy
git commit -s -m "Add executor.gateMaxWait config option"
```

---

## Task 3: Add monitor scaffolding — `readinessGates`, `gateExecutor`, `GateState`, `start()` init

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/processor/TaskPollingMonitor.groovy`
- Create: `modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy`

**Step 1: Write the failing test**

Create the new test file with the scaffolding test:

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * ...license header...
 */
package nextflow.processor

import java.util.concurrent.Future
import nextflow.Session
import nextflow.executor.ExecutorConfig
import spock.lang.Specification

class TaskPollingMonitorReadinessGateTest extends Specification {

    def 'should leave gateExecutor null when no gates are registered'() {
        given:
        def session = Mock(Session)
        def config = new ExecutorConfig([:])
        def monitor = new TaskPollingMonitor(name: 'local', session: session, config: config, pollInterval: '1s', capacity: 10)

        when:
        monitor.readinessGates = []   // simulate empty extension list

        then:
        monitor.gateExecutor == null
        monitor.gateStates.isEmpty()
    }

    def 'should expose GateState as a package-scoped type'() {
        when:
        def state = new TaskPollingMonitor.GateState([])
        then:
        state.futures == []
        state.scheduledAt > 0
    }
}
```

**Step 2: Run test to verify it fails**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: compilation error — fields and inner class do not exist.

**Step 3: Implement**

In `TaskPollingMonitor.groovy`:

Add imports near the top with the existing imports:

```groovy
import java.util.concurrent.Callable
import java.util.concurrent.ConcurrentHashMap
import java.util.concurrent.ConcurrentMap
import java.util.concurrent.ExecutionException
import java.util.concurrent.ExecutorService
import java.util.concurrent.Executors
import java.util.concurrent.Future
import groovy.transform.PackageScope
import nextflow.plugin.Plugins
import nextflow.util.CustomThreadFactory
import nextflow.util.Threads
```

Inside the class, alongside the existing fields:

```groovy
    @PackageScope
    List<TaskReadinessGate> readinessGates = Collections.emptyList()

    @PackageScope
    ExecutorService gateExecutor

    @PackageScope
    final ConcurrentMap<TaskHandler, GateState> gateStates = new ConcurrentHashMap<>()

    @PackageScope
    static class GateState {
        final long scheduledAt = System.currentTimeMillis()
        final List<Future<?>> futures
        GateState(List<Future<?>> futures) { this.futures = futures }
    }
```

In `start()`, add at the very top (before the existing barrier registration):

```groovy
        readinessGates = Plugins.getExtensions(TaskReadinessGate)
        if( readinessGates ) {
            gateExecutor = Threads.useVirtual()
                ? Executors.newVirtualThreadPerTaskExecutor()
                : Executors.newCachedThreadPool(new CustomThreadFactory('TaskReadinessGate'))
            session.onShutdown { gateExecutor.shutdownNow() }
            log.debug "Registered ${readinessGates.size()} task readiness gate(s): ${readinessGates*.class*.simpleName}"
        }
```

**Step 4: Run test to verify it passes**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: 2 tests pass.

Run also: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorTest'`
Expected: existing tests still pass.

**Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/TaskPollingMonitor.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy
git commit -s -m "Wire TaskReadinessGate discovery and executor into TaskPollingMonitor"
```

---

## Task 4: Submit gate work on `schedule()`

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/processor/TaskPollingMonitor.groovy:304-315`
- Modify: `modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy`

**Step 1: Write the failing test**

Add to `TaskPollingMonitorReadinessGateTest`:

```groovy
def 'should submit gate.prepare to executor on schedule'() {
    given:
    def latch = new CountDownLatch(1)
    def gate = Mock(TaskReadinessGate) {
        1 * prepare(_) >> { latch.countDown() }
    }
    def handler = mockHandler()
    def monitor = newMonitor()
    monitor.readinessGates = [gate]
    monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

    when:
    monitor.schedule(handler)

    then:
    latch.await(5, TimeUnit.SECONDS)
    monitor.gateStates.containsKey(handler)
    monitor.gateStates.get(handler).futures.size() == 1
}

def 'should not touch gate state when no gates are registered'() {
    given:
    def handler = mockHandler()
    def monitor = newMonitor()

    when:
    monitor.schedule(handler)

    then:
    monitor.gateStates.isEmpty()
}

private TaskPollingMonitor newMonitor() {
    new TaskPollingMonitor(name: 'local', session: Mock(Session), config: new ExecutorConfig([:]),
                           pollInterval: '1s', capacity: 10)
}

private TaskHandler mockHandler() {
    Mock(TaskHandler) { getTask() >> Mock(TaskRun) }
}
```

Add imports: `java.util.concurrent.CountDownLatch`, `java.util.concurrent.Executors`, `java.util.concurrent.TimeUnit`.

**Step 2: Run test to verify it fails**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: FAIL — `gateStates` empty because `schedule` does not populate it.

**Step 3: Implement**

Modify `schedule()` in `TaskPollingMonitor.groovy`:

```groovy
    @Override
    void schedule(TaskHandler handler) {
        pendingLock.lock()
        try{
            pendingQueue << handler
            if( readinessGates ) {
                final futures = new ArrayList<Future<?>>(readinessGates.size())
                for( gate in readinessGates ) {
                    futures << gateExecutor.submit({ gate.prepare(handler) } as Callable)
                }
                gateStates.put(handler, new GateState(futures))
            }
            taskAvail.signal()
            notifyTaskPending(handler)
            log.trace "Scheduled task > $handler"
        }
        finally {
            pendingLock.unlock()
        }
    }
```

**Step 4: Run test to verify it passes**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: all tests pass.

**Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/TaskPollingMonitor.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy
git commit -s -m "Submit TaskReadinessGate work on schedule()"
```

---

## Task 5: Wire `allGatesReady` into `canSubmit` — happy path

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/processor/TaskPollingMonitor.groovy:250-252` and add `allGatesReady`
- Modify: `modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy`

**Step 1: Write the failing test**

Add to `TaskPollingMonitorReadinessGateTest`:

```groovy
def 'should return false from canSubmit while gate prepare runs'() {
    given:
    def block = new CountDownLatch(1)
    def gate = Mock(TaskReadinessGate) {
        prepare(_) >> { block.await(5, TimeUnit.SECONDS) }
    }
    def handler = mockReadyHandler()
    def monitor = newMonitor()
    monitor.readinessGates = [gate]
    monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

    when:
    monitor.schedule(handler)
    def busy = monitor.canSubmit(handler)
    block.countDown()
    Thread.sleep(200)
    def ready = monitor.canSubmit(handler)

    then:
    !busy
    ready
    !monitor.gateStates.containsKey(handler)   // state removed on admission
}

private TaskHandler mockReadyHandler() {
    Mock(TaskHandler) {
        isReady() >> true
        canForkProcess() >> true
        getTask() >> Mock(TaskRun) { getName() >> 'task-x' }
    }
}
```

**Step 2: Run test to verify it fails**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: FAIL — `canSubmit` does not check gate readiness; returns true on first call.

**Step 3: Implement**

Modify `canSubmit` (line 250):

```groovy
    protected boolean canSubmit(TaskHandler handler) {
        allGatesReady(handler) \
            && handler.canForkProcess() \
            && handler.isReady() \
            && (capacity > 0 ? checkQueueCapacity(handler) : true)
    }

    private boolean allGatesReady(TaskHandler handler) {
        final state = gateStates.get(handler)
        if( !state ) return true

        for( f in state.futures ) {
            if( !f.isDone() ) return false
            if( f.isCancelled() ) {
                state.futures*.cancel(true)
                gateStates.remove(handler)
                throw new ProcessException("Task readiness gate was cancelled for task '${handler.task.name}'")
            }
            try { f.get() }
            catch( ExecutionException e ) {
                // cancel peer gates so their work doesn't outlive the failing task
                state.futures*.cancel(true)
                gateStates.remove(handler)
                throw e.cause instanceof ProcessException
                    ? (ProcessException) e.cause
                    : new ProcessException("Task readiness gate failed for task '${handler.task.name}'", e.cause)
            }
            catch( InterruptedException e ) {
                Thread.currentThread().interrupt()
                return false
            }
        }
        gateStates.remove(handler)
        return true
    }
```

Note: throwing from `canSubmit` is fine — the submit loop (`submitPendingTasks`) already wraps `canSubmit`/`submit` in a `try/catch(Throwable)` that calls `handleException(handler, e)` followed by `notifyTaskComplete(handler)`. Gate failures naturally route through the existing `errorStrategy` machinery. No new helper needed.

**Step 4: Run test to verify it passes**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: all tests pass.

Run also: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorTest'`
Expected: existing tests still pass.

**Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/TaskPollingMonitor.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy
git commit -s -m "Check TaskReadinessGate futures in canSubmit()"
```

---

## Task 6: Verify gate-exception routing (regression coverage)

> **Note:** This is a verification task. Task 5's `allGatesReady` already implements exception unwrapping; these tests pin that behavior in place. If they pass on the first run, that is the expected outcome — no implementation step.


**Files:**
- Modify: `modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy`

**Step 1: Write the failing test**

```groovy
def 'should throw ProcessException through canSubmit when gate throws'() {
    given:
    def gate = Mock(TaskReadinessGate) {
        prepare(_) >> { throw new ProcessException('boom') }
    }
    def handler = mockReadyHandler()
    def monitor = newMonitor()
    monitor.readinessGates = [gate]
    monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

    when:
    monitor.schedule(handler)
    awaitFutureCompletion(monitor, handler)
    monitor.canSubmit(handler)

    then:
    def e = thrown(ProcessException)
    e.message == 'boom'
    !monitor.gateStates.containsKey(handler)
}

def 'should wrap non-ProcessException causes'() {
    given:
    def gate = Mock(TaskReadinessGate) {
        prepare(_) >> { throw new RuntimeException('underlying') }
    }
    def handler = mockReadyHandler()
    def monitor = newMonitor()
    monitor.readinessGates = [gate]
    monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

    when:
    monitor.schedule(handler)
    awaitFutureCompletion(monitor, handler)
    monitor.canSubmit(handler)

    then:
    def e = thrown(ProcessException)
    e.cause.message == 'underlying'
}

private void awaitFutureCompletion(TaskPollingMonitor monitor, TaskHandler handler) {
    final state = monitor.gateStates.get(handler)
    for( f in state.futures ) {
        try { f.get(5, TimeUnit.SECONDS) }
        catch( Exception ignored ) { /* swallow — we only want to wait until done */ }
    }
}
```

**Step 2: Run tests**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: all tests pass. If they fail, the issue is in `allGatesReady`'s `ExecutionException` unwrap — fix the `e.cause instanceof ProcessException` branch.

**Step 3: Commit**

```bash
git add modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy
git commit -s -m "Test TaskReadinessGate exception propagation through canSubmit"
```

---

## Task 7: Enforce `gateMaxWait` timeout

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/processor/TaskPollingMonitor.groovy` (`allGatesReady`)
- Modify: `modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy`

**Step 1: Write the failing test**

```groovy
def 'should time out gate that exceeds gateMaxWait'() {
    given:
    def block = new CountDownLatch(1)
    def gate = Mock(TaskReadinessGate) {
        prepare(_) >> { block.await(10, TimeUnit.SECONDS) }
    }
    def handler = mockReadyHandler()
    def monitor = new TaskPollingMonitor(
        name: 'local', session: Mock(Session),
        config: new ExecutorConfig(gateMaxWait: '50 ms'),
        pollInterval: '1s', capacity: 10)
    monitor.readinessGates = [gate]
    monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()
    // gateMaxWait is wired from config inside the monitor constructor; no direct assignment needed.

    when:
    monitor.schedule(handler)
    Thread.sleep(150)
    monitor.canSubmit(handler)

    then:
    def e = thrown(ProcessException)
    e.message.contains('timed out')
    !monitor.gateStates.containsKey(handler)

    cleanup:
    block.countDown()
}
```

Plus a helper:

```groovy
def 'should not time out before gateMaxWait elapses'() {
    given:
    def block = new CountDownLatch(1)
    def gate = Mock(TaskReadinessGate) { prepare(_) >> { block.await() } }
    def handler = mockReadyHandler()
    def monitor = new TaskPollingMonitor(
        name: 'local', session: Mock(Session),
        config: new ExecutorConfig(gateMaxWait: '10 sec'),
        pollInterval: '1s', capacity: 10)
    monitor.readinessGates = [gate]
    monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

    when:
    monitor.schedule(handler)
    Thread.sleep(100)
    def busy = monitor.canSubmit(handler)
    block.countDown()
    Thread.sleep(100)
    def ready = monitor.canSubmit(handler)

    then:
    !busy
    ready
}
```

**Step 2: Run test to verify it fails**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: FAIL — timeout test does not see a ProcessException because the timeout check isn't implemented.

**Step 3: Implement**

Add a field to `TaskPollingMonitor` alongside the other gate fields:

```groovy
    @PackageScope
    Duration gateMaxWait
```

Initialise it in the constructor (around `TaskPollingMonitor.groovy:152-162`, after `this.config = params.config as ExecutorConfig`):

```groovy
    this.gateMaxWait = config?.gateMaxWait
```

The `config?` guard mirrors the defensive style of the surrounding lines (some tests construct the monitor without a config). For the production path, `config.gateMaxWait` is non-null because Task 2 wired a 24h default.

In `allGatesReady`, add a timeout check before the future loop:

```groovy
        if( gateMaxWait && System.currentTimeMillis() - state.scheduledAt > gateMaxWait.toMillis() ) {
            state.futures*.cancel(true)
            gateStates.remove(handler)
            throw new ProcessException("Task readiness gate timed out after ${gateMaxWait} for task '${handler.task.name}'")
        }
```

Place this check immediately after the `if( !state ) return true` early-return, before the future-iteration loop.

**Step 4: Run test to verify it passes**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: all tests pass.

**Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/TaskPollingMonitor.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy
git commit -s -m "Enforce executor.gateMaxWait timeout in TaskReadinessGate flow"
```

---

## Task 8: Cancel gate futures on `evict()`

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/processor/TaskPollingMonitor.groovy:327-344`
- Modify: `modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy`

**Step 1: Write the failing test**

```groovy
def 'should cancel gate futures on evict'() {
    given:
    def started = new CountDownLatch(1)
    def interrupted = new CountDownLatch(1)
    def gate = new TaskReadinessGate() {
        void prepare(TaskHandler h) throws InterruptedException {
            started.countDown()
            try { new CountDownLatch(1).await() }
            catch( InterruptedException e ) { interrupted.countDown(); throw e }
        }
    }
    def handler = mockReadyHandler()
    def monitor = newMonitor()
    monitor.readinessGates = [gate]
    monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

    when:
    monitor.schedule(handler)
    started.await(5, TimeUnit.SECONDS)
    monitor.evict(handler)

    then:
    interrupted.await(5, TimeUnit.SECONDS)
    !monitor.gateStates.containsKey(handler)
}
```

**Step 2: Run test to verify it fails**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: FAIL — `interrupted` latch never trips because `evict` doesn't cancel.

**Step 3: Implement**

Modify `evict()`:

```groovy
    @Override
    boolean evict(TaskHandler handler) {
        if( !handler ) {
            return false
        }

        gateStates.remove(handler)?.futures*.cancel(true)

        if( remove(handler) ) {
            // ... existing logic unchanged
        }
        return false
    }
```

**Step 4: Run test to verify it passes**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: all tests pass.

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorTest'`
Expected: existing tests still pass.

**Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/TaskPollingMonitor.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy
git commit -s -m "Cancel in-flight TaskReadinessGate futures on task eviction"
```

---

## Task 9: Verify multi-gate coverage (regression tests)

> **Note:** Verification task. Tasks 4-7 already implement the multi-gate semantics; these tests pin them.


**Files:**
- Modify: `modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy`

**Step 1: Write the failing test**

```groovy
def 'should require all gates to complete before admission'() {
    given:
    def latch1 = new CountDownLatch(1)
    def latch2 = new CountDownLatch(1)
    def g1 = Mock(TaskReadinessGate) { prepare(_) >> { latch1.await() } }
    def g2 = Mock(TaskReadinessGate) { prepare(_) >> { latch2.await() } }
    def handler = mockReadyHandler()
    def monitor = newMonitor()
    monitor.readinessGates = [g1, g2]
    monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

    when:
    monitor.schedule(handler)
    Thread.sleep(50)
    def busy0 = monitor.canSubmit(handler)
    latch1.countDown()
    Thread.sleep(50)
    def busy1 = monitor.canSubmit(handler)
    latch2.countDown()
    Thread.sleep(50)
    def ready = monitor.canSubmit(handler)

    then:
    !busy0
    !busy1
    ready
}

def 'should fail fast when one of several gates throws'() {
    given:
    def block = new CountDownLatch(1)
    def g1 = Mock(TaskReadinessGate) { prepare(_) >> { throw new ProcessException('g1 failed') } }
    def g2 = Mock(TaskReadinessGate) { prepare(_) >> { block.await() } }
    def handler = mockReadyHandler()
    def monitor = newMonitor()
    monitor.readinessGates = [g1, g2]
    monitor.gateExecutor = Executors.newVirtualThreadPerTaskExecutor()

    when:
    monitor.schedule(handler)
    Thread.sleep(100)
    monitor.canSubmit(handler)

    then:
    thrown(ProcessException)
    // g2 is left in flight; the eviction path in the production submit loop will cancel.

    cleanup:
    block.countDown()
}
```

**Step 2: Run tests**

Run: `./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'`
Expected: both tests pass. The "fail fast" test relies on `isDone()` returning true for the failed gate before the still-running one — true under all reasonable scheduling because the thrown gate completes immediately. If it flakes, add `awaitFutureCompletion(monitor, handler)` before the `canSubmit` call.

**Step 3: Commit**

```bash
git add modules/nextflow/src/test/groovy/nextflow/processor/TaskPollingMonitorReadinessGateTest.groovy
git commit -s -m "Test multi-gate semantics for TaskReadinessGate"
```

---

## Task 10: Run the full monitor test class to catch regressions

**Step 1: Run**

```bash
./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorTest'
./gradlew :nextflow:test --tests 'nextflow.processor.TaskPollingMonitorReadinessGateTest'
./gradlew :nextflow:test --tests 'nextflow.executor.ExecutorConfigTest'
```

Expected: all green.

**Step 2: Run the broader `nextflow.processor` test package**

```bash
./gradlew :nextflow:test --tests 'nextflow.processor.*'
```

Expected: all green. Investigate any new failures — likely indicates a behavioural regression introduced by the `canSubmit` reordering or the new fields.

**Step 3: Commit (only if any incidental cleanup needed)**

If tests reveal nothing, no commit. Otherwise fix and commit with a descriptive message.

---

## Task 11: Documentation

**Files:**
- Create: `docs/developer/task-readiness-gate.md` (or append to an existing developer-docs page if conventions dictate — locate first via `ls docs/developer/`)
- Modify: `docs/reference/config.md` (find the `executor` scope section and add the `gateMaxWait` entry)

**Step 1: Write the developer page**

Content roughly:

```markdown
(task-readiness-gate)=

# TaskReadinessGate

`TaskReadinessGate` is a plugin extension point that defers task submission until an external precondition is met — for example, restoring an S3 object from Glacier before an AWS Batch worker tries to stage it.

## Contract

A gate implements one method:

```groovy
package nextflow.processor

interface TaskReadinessGate extends ExtensionPoint {
    void prepare(TaskHandler handler) throws InterruptedException
}
```

- **Blocking is allowed.** `prepare` runs on a managed virtual-thread executor; calling `Thread.sleep`, blocking I/O, or long-polling APIs is fine. The scheduler thread is never blocked.
- **Throwing fails the task.** Any exception marks the task as permanently failed with the thrown cause routed through the `errorStrategy` directive. `ProcessException` is canonical.
- **Interrupts must be honored.** Task eviction, workflow abort, and the `executor.gateMaxWait` backstop cancel the in-flight `prepare` by interrupting its thread. Use interruptible primitives.
- **Multiple gates compose.** When several plugins register gates, all must complete successfully before the task is admitted. Gates run in parallel; order is unspecified.

## Example

```groovy
class GlacierReadinessGate implements TaskReadinessGate {
    @Override
    void prepare(TaskHandler handler) throws InterruptedException {
        if( handler.task.config.hints['glacier/skip'] == true ) return
        // issue RestoreObject, poll HeadObject until restored, throw on permanent failure
    }
}
```

## Configuration

`executor.gateMaxWait` bounds the time a gate may spend preparing a single task (default `24h`). Set to `null` to wait indefinitely.

## Per-process opt-out

Use the `hints` directive:

```nextflow
process FAST_PATH {
    hints 'glacier/skip': true
    // ...
}
```

Gates inspect `handler.task.config.hints` and return early when set.
```

**Step 2: Add config reference entry**

Locate the `executor.*` section in `docs/reference/config.md` and add `gateMaxWait` alphabetically.

**Step 3: Commit**

```bash
git add docs/
git commit -s -m "Document TaskReadinessGate extension point and gateMaxWait config"
```

---

## Task 12: Changelog

**Files:**
- Modify: `changelog.txt`

**Step 1: Add entries**

Open `changelog.txt`. Entries land as flat bullet lines under the unreleased version header (no `NEW FEATURES` subsection — see existing format). Add two lines under the topmost version block:

```
- Add TaskReadinessGate plugin extension point for deferring task submission until external preconditions are met
- Add executor.gateMaxWait config option as a safety net for stuck readiness gates
```

The `[hash]` suffix that appears on existing lines is added by the release process; do not add one manually.

**Step 2: Commit**

```bash
git add changelog.txt
git commit -s -m "Add changelog entries for TaskReadinessGate and executor.gateMaxWait"
```

---

## Task 13: Final verification

**Step 1: Compile cleanly**

```bash
./gradlew :nextflow:compileGroovy
```

Expected: BUILD SUCCESSFUL.

**Step 2: Full nextflow-module test run**

```bash
./gradlew :nextflow:test
```

Expected: BUILD SUCCESSFUL. This is longer (~5-10 min) but verifies no regressions in unrelated tests.

**Step 3: Push and open PR**

```bash
git push -u origin feature/task-readiness-gate
gh pr create --base master --title "Add TaskReadinessGate plugin extension point" --body "$(cat <<'EOF'
## Summary

Implements the `TaskReadinessGate` plugin extension point described in ADR `20260516-task-readiness-gate.md` (#7151), using the blocking-`prepare` design that emerged from review of that ADR.

Plugins implement `void prepare(TaskHandler) throws InterruptedException`. The method may block freely; core runs it on a managed virtual-thread executor and polls the resulting future from `canSubmit`. Throwing routes through the existing `errorStrategy` machinery. An `executor.gateMaxWait` option (default `24h`) bounds wait time as a safety net for stuck gates.

Behavior is bit-identical when no plugin registers a gate.

## Test plan

- [ ] CI: full `:nextflow:test` green
- [ ] Spock specs cover: no-gates path, happy-path single gate, gate throws, gate timeout, multi-gate AND semantics, eviction cancellation
- [ ] `TaskPollingMonitorTest` (existing) still green
EOF
)"
```

Expected: PR URL returned.

---
