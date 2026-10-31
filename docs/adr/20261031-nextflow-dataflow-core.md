# Replace GPars with a Nextflow-owned dataflow core

- Authors: Ben Sherman
- Status: accepted
- Deciders: Ben Sherman
- Date: 2026-10-31
- Tags: dataflow, runtime, performance, dependencies

## Summary

Absorb the part of GPars that Nextflow uses into a new `nf-dataflow` module, keep its package and class names, and replace its internals: serial operators on a shared pool instead of actors, lightweight channels, and backpressure at processes and sources.

## Problem Statement

Nextflow builds every channel, operator and process on GPars 1.2.1, which has had no release since 2014. Around 90 Nextflow files and several third-party plugins import `groovyx.gpars`.

- **Memory.** A buffered channel item costs about 290 B (a `DataflowVariable`, a stream node, an event orchestrator and their lists). Nothing bounds the number of buffered items or pending tasks: a fast source feeding a slow process runs out of memory before the first task completes, and 10⁶ items through five operators do not fit in a 256 MiB heap.
- **Overhead.** Operators are actors: each firing goes through a mailbox message, a `HashMap` and a sort, and re-registers a `whenBound` callback per input.
- **Maintenance.** The actor framework, `jsr166y` and `multiverse` are dead weight that Nextflow cannot fix upstream.

## Goals or Decision Drivers

- Bound the memory of a run when tasks are produced faster than they are executed.
- Keep the semantics that scripts, Nextflow and plugins rely on: per-operator serial bodies, `maxForks`, poison pills, event listeners, blocking `getVal()`.
- Keep the classes that Nextflow and its core plugins use; other plugins may need a rebuild.
- Keep blocking work on dataflow threads safe: task preparation (staging, hashing, cache DB), I/O operators (`split*`, `collectFile`) and user closures.
- Small, incremental changes on the Nextflow side.

## Non-goals

- A new dataflow API or package name.
- Splitting task execution from the process operator.

## Considered Options

- Keep GPars as is.
- A single dataflow loop that runs all operators on one thread.
- Serial operators on a shared pool, in a trimmed fork of GPars.
- An external reactive library (Reactor) or virtual threads with bounded queues.

## Pros and Cons of the Options

### Keep GPars

- Good, because nothing changes.
- Bad, because memory stays unbounded and the overhead stays per item.
- Bad, because backpressure would have to be retrofitted onto unmaintained actor code anyway.

### Single dataflow loop

- Good, because it is the fastest option and deterministic.
- Bad, because user closures can do blocking I/O (e.g. `file(..., checkIfExists: true)` on remote storage), which would stall every operator and all task creation. Task preparation and the built-in I/O operators could be moved off the loop, but arbitrary closures cannot.

### Serial operators on a shared pool

- Good, because each operator keeps GPars' guarantee of one body at a time, so operator state stays safe without changes.
- Good, because a blocking closure only delays its own operator, as today.
- Good, because a pool of one thread is the single loop, which gives a deterministic test mode.
- Bad, because handing items between operators goes through the pool's task queue (amortized over batches of up to 128 items), which is slower than direct calls on a single loop.

### Reactor or bounded queues

- Bad, because a bounded channel or demand-based backpressure deadlocks a process that reads a channel and a value reduced from it (`P(ch, ch.collect())`).
- Bad, because every operator would have to be rewritten against a new model, and Reactor adds a dependency.

## Solution or decision outcome

Vendor GPars as `modules/nf-dataflow` (Groovy as its only dependency), trimmed to what Nextflow and its plugins use, and replace its internals behind the same names. Operators run as mailboxes on a shared pool. Backpressure is on by default.

## Rationale & discussion

**Absorbed and trimmed.** The fork keeps the `groovyx.gpars` package and the names Nextflow and its plugins use, and trims or replaces everything else:

| GPars API | Outcome |
|---|---|
| `DataflowVariable`, `DataflowQueue`, `DataflowBroadcast`, `DataflowReadChannel`, `DataflowWriteChannel`, `DataflowProcessor`, `PoisonPill`, `Agent` | Kept (rewritten internally), so most importing files do not change |
| Actor framework | Replaced by the mailbox engine |
| `Dataflow` facade, parallel groups (`PGroup`) | Replaced by the shared pool |
| `DataflowExpression` | Merged into `DataflowVariable` |
| `DataflowEventAdapter` | Merged into `DataflowEventListener` (default methods) |
| Broadcast writer adapter | Merged into `DataflowBroadcast` |
| Broadcast reader (`stream.DataflowStreamReadAdapter`) | Moved to `groovyx.gpars.dataflow` |
| `leftShift` (`<<`) | Replaced by `bind` (`send` for `Agent`) |
| `ControlMessage` | Replaced by `PoisonPill` (the type of `Channel.STOP`) |
| Derived values (unknown property or method of a value) | Removed: fails with an error suggesting an operator |
| `Promise`, timed reads, `DataflowChannel`, unused processor methods and listener events | Removed |
| Remote and serial actors, STM, CSP, `GParsPool`, selectors, `jsr166y`, `multiverse` | Removed |

The renames changed a few Nextflow call sites mechanically. GPars types and methods are not script API, so most of these changes only affect plugins that use them. The exception is derived values: a script that calls a non-operator method or property on a value channel (e.g. `.toList().sort()`) now fails. A JUnit conformance suite pins the GPars semantics that Nextflow relies on. The steps landed as separate commits: vendor with no behavior change, backpressure behind a flag, the mailbox engine with the actor engine as a fallback, the channel rewrite, removal of the fallback, then the API reduction.

**Concurrency model.** Each `DataflowProcessor` is a mailbox. When an input hands over a value, the operator is scheduled on the shared pool unless it is already scheduled, and the drain fires on complete input sets in a bounded batch before yielding the thread. A read takes a buffered value at once or subscribes for the next one: no message or `DataflowVariable` per item. With `maxForks > 1` the body is forked onto the pool, up to `maxForks` at once; when all are busy, the operator stops reading inputs instead of blocking a thread. `NXF_DATAFLOW_POOL_SIZE` fixes the pool size. With `1`, the network runs on one loop thread in a deterministic order, and a blocking read on that thread fails instead of hanging.

**Channels.** A queue is a deque with a list of waiting readers. A broadcast is a linked list of nodes shared by its readers (24 B per item). A value is a write-once cell. `getVal()` still blocks, so tests and blocking reads keep working.

**Backpressure.** Operators and channels never block; the bound sits at processes and sources:

- A process stops taking new inputs while the tasks it is preparing plus its tasks waiting to be submitted reach the high-water mark (`NXF_DATAFLOW_HIGH_WATER`, default `max(64, 2 × CPUs)`), and resumes at half of it.
- A source factory (`of`, `fromList`, `fromPath`, `fromFilePairs`, `interval`, `watchPath`) blocks its own producer thread while any process downstream of it is over the mark.

Since nothing on a dataflow thread waits for a process, every pause drains, and items waiting for an unbound value input are not counted, so `P(ch, ch.collect())` cannot deadlock. An abort or cancel opens every gate. There is no runtime safety net: a lost backlog decrement would show up as a hang, so the tests check that every run ends with each process backlog back to zero and no source blocked, including on retried, ignored, skipped, cached and failed tasks. `NXF_DATAFLOW_HIGH_WATER=0` disables backpressure and restores the unbounded behavior; a negative or non-integer value fails the run at startup.

**Results.** On the benchmarks of `modules/nf-dataflow/src/bench`: 10⁶ items through five operators run in a 256 MiB heap (GPars: out of memory), with a third of the peak RSS; a fast source feeding a slow process completes in 128 MiB with a flat live heap (GPars: out of memory before the first task). Process-heavy runs are 10–15% faster.

## Consequences

Nextflow has no migration page for the next release yet. These notes should go in it:

- Calling a non-operator method or property on a value channel, such as `ch.toList().sort()`, now fails. Use `map` instead, e.g. `ch.toList().map { it.sort() }`.
- Plugins that use the GPars APIs `leftShift` (`<<`), `DataflowEventAdapter`, `ControlMessage`, `DataflowExpression`, `PGroup` or `Dataflow`, or the `groovyx.gpars.dataflow.stream` package, must be updated and rebuilt.
- Backpressure is enabled by default. Set `NXF_DATAFLOW_HIGH_WATER=0` to restore the previous unbounded behavior.

## Links

- The single-loop engine and the rule of bounding processes and sources, never operator queues, come from an internal Seqera prototype.
- [GPars](https://github.com/GPars/GPars), the vendored upstream (1.2.1).
