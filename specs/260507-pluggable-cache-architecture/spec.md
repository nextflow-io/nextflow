# Feature Specification: Pluggable cache architecture for global cache implementations

**Feature Branch**: `20260202-global-cache` (continuation)
**Created**: 2026-05-07
**Status**: Draft
**Input**: Refactor Nextflow's task caching to define plugin-extensible interfaces that allow different global cache implementations to be provided as plugins. Today's behaviour stays as the default implementation. Builds on PR #6927 (versioned task hasher) and `adr/20260202-global-cache.md`.

## Motivation

Today's task caching is hardcoded across `TaskProcessor`, `TaskHasher`, `CacheDB`, `PublishDir`, and `FilePorter`. The workdir doubles as both execution scratch space and durable output store. The hasher is a concrete class with hardcoded keys. The cache lookup is a tight loop in `TaskProcessor` that interleaves `LockManager`-based session-local serialization, hash bumping for collision avoidance, and physical workdir scans for `.exitcode`, `.command.out`, and the output files via `collectOutputs`.

This entanglement makes it impossible to deliver alternative cache implementations — including the global cache described in `adr/20260202-global-cache.md` — without forking core code. Different cache implementations need different decisions about:

- **Task identity**: which keys go into the hash; how files are fingerprinted (path+attrs vs. deep content vs. cached deep hash).
- **Coordination**: who owns the slot for a given hash; whether collisions wait, bump, or fail.
- **Storage**: where output files physically live after a task completes; whether the workdir survives.
- **Reference counting**: which inputs and outputs are still in use; when files are safe to delete.
- **Cleanup policy**: TTL, ref-counted, manual; in-process or out-of-band.

A plugin-extensible architecture moves each of these decisions behind an interface that cache implementations override, while the default implementation preserves today's behavior byte-for-byte.

## Goals

- Extract five plugin-extensible seams in cache-related logic: task hashing, cache resolution (lookup + claim), output adoption, file-usage events, optional cleanup.
- Move all cache-related decisions out of `TaskProcessor` and into the cache interface.
- Enable the global cache (and any future cache variant) to ship as a plugin without forking core code.
- Preserve byte-identical behavior when no plugin is registered and no opt-in config is set.
- Keep on-disk format of LevelDB and `nf-cloudcache` unchanged.

## Non-Goals

- No change to executor APIs (`TaskHandler`, `Executor` subclasses).
- No change to `task.getWorkDirFor(hash)` or workdir naming conventions.
- No replacement of `TraceObserverV2`. Cache events are independent.
- No change to `storeDir` semantics; `storeDir`-based outputs continue to bypass the cache entirely.
- No per-process hasher selection (single hasher per session).
- No async or batched event emission. All events are synchronous.
- No change to error/retry strategies (`maxErrors`, `maxRetries`, dynamic resources).

## Architecture Overview

```
                     ┌──────────────────────────┐
                     │     TaskProcessor        │
                     │  (orchestration only)    │
                     └────────────┬─────────────┘
                                  │
            ┌─────────────────────┼─────────────────────┐
            │                     │                     │
            ▼                     ▼                     ▼
   ┌────────────────┐    ┌────────────────┐    ┌──────────────────┐
   │  TaskHasher    │    │   Cache        │    │ FilePorter /     │
   │  (interface,   │    │ (interface ext)│    │  PublishDir      │
   │   pluggable)   │    │                │    │ (call notify*)   │
   └───────┬────────┘    └───────┬────────┘    └─────────┬────────┘
           │                     │                       │
           └──────────┬──────────┴───────────────────────┘
                      ▼
        ┌──────────────────────────────────────┐
        │       Cache Plugin Implementation    │
        │  - TaskHasher                        │
        │  - CacheStore (KV state)             │
        │  - TaskCacheEntry restore            │
        │    (outputs-shaped)                  │
        │  - adopt(task, outputs)              │
        │  - beginTask / endTask events        │
        │  - notifyPublish / notifyFilePort    │
        │  - (optional) CacheReader            │
        │  - (optional) CacheCleaner           │
        │  - (optional) CacheDropper           │
        └──────────────────────────────────────┘
```

The five seams:

1. **Pluggable `TaskHasher`** — interface + factory (extending PR #6927). Plugin-registered hashers; cache implementation declares its preferred default.
2. **Outputs-shaped cache contract** — `Cache.getTaskCacheEntry()` returns resolved output bindings (URIs + values + exit code + stdout); `TaskProcessor` no longer scans the workdir for restore. Default implementation realizes bindings lazily via a workdir scan at restore time.
3. **Workdir adoption hook** — `cache.adopt(task, outputs)` after task success; default no-op; cache returns `WorkdirDisposition` (KEEP/DELETE). Workdir derivation (`task.getWorkDirFor(hash)`) untouched.
4. **File-usage event methods on the cache** — task lifecycle (`beginTask`/`endTask`), publish, and file porting events extending the cache interface; default no-ops; events are authoritative source of truth for ref counts.
5. **Capability interfaces for read/cleanup outside `Cache`** — record iteration is exposed via an optional `CacheReader`; selective bulk eviction via `CacheCleaner`; unconditional destructive ops via `CacheDropper`; the composite `CacheCleanup` is the return type of `Cache.openForClean()`. `nextflow log` / `nextflow cache clean` are thin pass-throughs that select the right capability via `openForRead()` / `openForClean()` and refuse gracefully when the implementation does not support a verb.

## Detailed Design

### Section 1 — Pluggable `TaskHasher`

Builds on PR #6927 which already extracted `TaskHasher` to an interface with `AbstractTaskHasher` + `TaskHasherV1`/`V2` and a `Version`-keyed factory cached on the `Session`. Selection is extended from "env-var-driven enum" to "plugin-driven name resolution"; the cache implementation drives the default.

**Interface (from #6927, unchanged):**

```groovy
interface TaskHasher {
    HashCode compute()
    Map<String,Object> getTaskGlobalVars()
    List<Path> getTaskBinEntries(String script)
}

abstract class AbstractTaskHasher implements TaskHasher { /* shared helpers */ }
```

**Factory (replaces `TaskHasherFactory.Version` with a name registry):**

```groovy
interface TaskHasherFactory {
    /** Stable identifier ("v1", "v2", "global", etc.) */
    String getName()
    TaskHasher create(TaskRun task)
}
```

Plugins register `TaskHasherFactory` instances via the standard plugin extension mechanism, mirroring how `CacheFactory` is plugin-discovered today. Built-ins: `v1`, `v2`.

**Resolution precedence** (evaluated once at session init, cached on `Session.hashStrategy`):

1. Cache-driven default: the active `CacheStore` implementation declares a preferred hasher via `CacheFactory.getDefaultTaskHasher()` (returns name).
2. Session-level config: `nextflow.config` top-level `taskHasher` setting.
3. Env var: `NXF_TASK_HASH_VER` (kept for backward-compat with #6927).
4. Hardcoded default: `v2`.

File-hash strategy stays internal to the hasher implementation. A hypothetical `GlobalTaskHasher` calls `CacheHelper.hasher(value, HashMode.DEEP, version)` internally; it is not a third pluggable axis. Backward compat: with no plugin and no config, resolution lands on `v2` and produces byte-identical hashes to today.

### Section 2 — Cache Resolution: `beginTask` / `endTask`

The cache becomes the single owner of all per-task lifecycle decisions: hit lookup, slot claim, collision resolution, and input ref-counting. `TaskProcessor` makes one call and gets back one of two outcomes.

```groovy
interface Cache extends Closeable {

    /** Initialise for active (writable) use. Returns self for chaining. */
    Cache open()

    /**
     * Initialise the underlying store in read-only mode and narrow the returned
     * handle to the read interface. The runtime object MAY be the same instance
     * (with internal state tracking the mode). Throws {@link UnsupportedOperationException}
     * when the implementation does not support browsing.
     */
    CacheReader openForRead()

    /**
     * Initialise for cleanup operations and narrow to the {@link CacheCleanup}
     * composite. Throws {@link UnsupportedOperationException} when the
     * implementation does not support cleanup.
     */
    CacheCleanup openForClean()

    // close() inherited from Closeable

    /**
     * Resolve the cache state for a task about to start. Single entry point for:
     *   - cache-hit lookup
     *   - collision resolution + claim acquisition
     *   - input ref-counting (cache impl uses `inputs` for its ledger)
     */
    CacheResolution beginTask(BeginTaskRequest req)

    /**
     * Called once after task lifecycle ends (success or failure). Single entry point for:
     *   - claim release
     *   - input ref-decrement
     *   - retry-friendly state cleanup on failure
     */
    void endTask(EndTaskRequest req)

    AdoptionResult adopt(TaskRun task, Map<String, OutputBinding> resolvedOutputs)
    void putTaskCacheEntry(TaskHandler handler, TraceRecord trace, Map<String, OutputBinding> outputs)
    TaskCacheEntry getTaskCacheEntry(HashCode hash, TaskProcessor processor)

    void notifyPublish(URI source, URI destination, PublishMode mode)
    void notifyFilePort(URI remote, URI local, TaskRun task)
}

class BeginTaskRequest {
    HashCode hash             // unbumped hash from TaskHasher
    TaskRun task
    int retryCount            // task.failCount
    List<URI> inputs          // post-port input file URIs (for ref counting)
    boolean tryHit            // false for retry paths that should skip hit lookup
}

class EndTaskRequest {
    HashCode finalHash        // hash actually claimed (from beginTask result)
    TaskRun task
    boolean success
}

class CacheResolution {
    enum Outcome { CACHE_HIT, CLAIM_GRANTED }
    Outcome outcome
    HashCode finalHash
    TaskCacheEntry entry      // populated when outcome == CACHE_HIT
}
```

The cache implementation owns the entire resolve-or-claim loop, including the hash-bump policy used for collision avoidance. Different cache implementations may adopt different policies (today's bump-and-retry, waiting, strict de-duplication, exponential backoff) without `TaskProcessor` changes.

**Updated `TaskProcessor.checkCachedOrLaunchTask`:**

```groovy
final protected void checkCachedOrLaunchTask(TaskRun task, HashCode hash, boolean shouldTryCache) {
    final inputs = collectInputUris(task)
    final result = session.cache.beginTask(new BeginTaskRequest(
        hash:       hash,
        task:       task,
        retryCount: task.failCount,
        inputs:     inputs,
        tryHit:     shouldTryCache
    ))

    if (result.outcome == CACHE_HIT) {
        restoreFromCacheEntry(task, result.entry, result.finalHash)
        return
    }
    submitTask(task, result.finalHash, task.getWorkDirFor(result.finalHash))
}
```

`TaskProcessor.finalizeTask` unconditionally calls `session.cache.endTask(...)` after the existing finalization.

**Default implementation** (`DefaultCache`) preserves today's behavior verbatim. It wraps a `CacheDB` instance (the existing class, untouched) and a `LockManager`:

- `beginTask`: implements the existing `while` loop — bump hash by `tries`, attempt to load `TaskEntry` via the wrapped `CacheDB`, check `LockManager.acquire` + `workDir.exists()` + `mkdirs()`. Returns `CACHE_HIT` or `CLAIM_GRANTED`.
- `endTask`: no-op for ref counting (the local lock was released inside `beginTask`).
- `getTaskCacheEntry`: delegates to the wrapped `CacheDB.getTaskEntry`, opens the workdir, scans for `.exitcode`/`.command.out`/output files, assembles `TaskCacheEntry`.
- `putTaskCacheEntry`: delegates to the wrapped `CacheDB`'s existing async writer.
- `adopt`: returns input bindings unchanged with `WorkdirDisposition.KEEP`.

`LockManager` moves out of `TaskProcessor` and becomes an internal of the `DefaultCache` implementation. This is the only meaningful relocation.

### Section 3 — Outputs-Shaped Cache Restore

Today the workdir is consulted four times during cache restore: `getTaskEntry` (DB), `resumeDir.exists()` (workdir), `.exitcode` (workdir), `collectOutputs` (workdir). The outputs-shaped contract replaces those with a single typed lookup that returns everything needed to bind output channels.

```groovy
class TaskCacheEntry {
    HashCode hash
    TraceRecord trace
    TaskContext context              // existing — needed for cacheable values
    Integer exitCode
    URI stdoutUri                    // points at stored stdout
    Map<String, OutputBinding> outputs    // keyed by output param name
}

class OutputBinding {
    OutputType type                  // FILE | VALUE | EVAL
    Object value                     // primitive value, or List<URI> for files, or eval string
    Map<String, Object> metadata     // size, mtime, content hash — opaque to TaskProcessor
}
```

**Default implementation behavior:**

- `putTaskCacheEntry`: stores `TaskEntry` exactly as today (trace + context + ref count). Output bindings are NOT serialized. On-disk format is unchanged.
- `getTaskCacheEntry`: loads the `TaskEntry`, opens the workdir from `trace.workDir`, reads `.exitcode`, `.command.out`, and runs `collectOutputs(task)` to assemble the `OutputBinding` map. The workdir scan moves *inside* the cache implementation rather than living in `TaskProcessor`.

**Alternative implementations** serialize output bindings into their own store on `putTaskCacheEntry` and return them directly on `getTaskCacheEntry`, with no workdir scan.

**`TaskProcessor` consequence:** `collectOutputs(task)` and the manual `.exitcode`/`.command.out` reads are removed from the cache-restore path. The new `restoreFromCacheEntry(task, entry, hash)` helper consumes a `TaskCacheEntry` and binds outputs directly. The success path (`finalizeTask`) and `checkStoredOutput`-for-`storeDir` keep calling `collectOutputs` as today.

### Section 4 — Workdir Adoption

After successful task execution, the cache may adopt outputs into its own store and signal whether the workdir is disposable.

```groovy
class AdoptionResult {
    Map<String, OutputBinding> outputs    // possibly rewritten URIs
    WorkdirDisposition disposition        // KEEP | DELETE
}

enum WorkdirDisposition { KEEP, DELETE }
```

**Call ordering inside `TaskProcessor.finalizeTask`** (success path):

```
1. collectOutputs(task)                                          // resolves OutputBinding from workdir
2. AdoptionResult adopted = cache.adopt(task, outputs)           // NEW
3. cache.putTaskCacheEntry(handler, trace, adopted.outputs)
4. bind adopted.outputs to output channels                       // downstream sees rewritten URIs
5. cache.endTask(EndTaskRequest{ success: true, ... })
6. if (adopted.disposition == DELETE) deleteWorkdirAsync(task.workDir)
```

Adoption happens *before* channel binding so downstream tasks consume the rewritten URIs, and *before* persisting the cache entry so future hits return canonical (cache-managed) URIs. Workdir deletion is async and last; failures during delete are logged but do not fail the task.

**Contract — DELETE invariant:** if a cache implementation returns `WorkdirDisposition.DELETE`, all FILE-typed bindings MUST point outside the workdir AND `TaskCacheEntry.stdoutUri` (when later retrieved) MUST resolve outside the workdir. The cache implementation is responsible for replicating stdout into managed storage before returning `DELETE`. Both are enforced by assertions at adopt time; violations fail the run loudly.

**Default implementation:** returns input bindings unchanged with `KEEP`. Behavior identical to today.

**Failed tasks do not adopt.** `adopt()` is only called on success. On failure, `endTask(success=false)` fires; failure-path workdir handling follows today's existing error-strategy semantics. The cache implementation can do its own cleanup of partial state via the failure `endTask` event.

### Section 5 — File-Usage Tracking

Cache implementations that perform reference-counted cleanup require visibility into file consumer/producer relationships beyond task lifecycle. Two standalone events cover the remaining cases.

| Method | Call site | When |
|--------|-----------|------|
| `notifyPublish` | `PublishDir.apply` / `PublishOp` | After successful publish operation |
| `notifyFilePort` | `FilePorter.transfer` | After successful download |

Default implementations are empty methods. LevelDB and `nf-cloudcache` continue to perform no ref counting; behavior is unchanged.

**Authoritative-via-events contract** — the things this design commits to:

1. **`beginTask` / `endTask` are paired exactly once per task lifecycle.** Every code path that fires `beginTask` MUST guarantee `endTask` fires, including: successful execution, task failure (any error strategy), task termination via session abort, dynamic resources retry (counts as one start/one complete; the retry is a new task with its own pair), JVM shutdown (best-effort via shutdown hook, documented as a known gap).
2. **`endTask` fires AFTER channel binding completes.** This guarantees that when downstream tasks get their `beginTask`, the upstream task's "still in use" ref count reflects the new consumer before the upstream attempts to release. Without this ordering, a cache implementation could see ref=0 momentarily between upstream completing and downstream starting.
3. **`inputs` in `BeginTaskRequest` are post-port URIs** — what the task will physically read.
4. **Cache-hit ref counting is internal.** When `beginTask` returns `CACHE_HIT`, the cache implementation produces that return value; it knows the hit happened and updates ref counts inline before returning. No separate `notifyCacheHit` event.
5. **No retroactive corrections.** Events represent reality at the moment of emission. Cache state is built incrementally; there is no "rebuild from history" path.

### Section 6 — Optional capability interfaces

Beyond the mandatory `Cache` lifecycle, alternative caches opt into three orthogonal capabilities. Each is a separate Groovy interface in `nextflow.cache`; a plugin declares `implements` for exactly the ones it supports. The CLI dispatches via `instanceof` checks after asking `Cache.openForRead()` / `openForClean()` to narrow the handle's type.

#### `CacheReader` — record read + iteration

```groovy
interface CacheReader extends Closeable {
    TraceRecord getTraceRecord(HashCode hash)
    TraceRecord findTraceRecord(Closure<Boolean> criteria)
    void eachRecord(Closure consumer)
}
```

Used by `nextflow log <runName>` and any per-run reporting. The receiver is obtained via `cache.openForRead()`, which returns the (narrowed) reader handle; the underlying store opens in read-only mode (no claim ownership, no async writers).

#### `CacheCleaner` — selective bulk cleanup

```groovy
interface CacheCleaner {
    CleanupReport cleanAll(boolean dryRun)
    CleanupReport cleanByHash(HashCode hash, boolean dryRun)
    CleanupReport cleanOlderThan(Duration age, boolean dryRun)
    CleanupReport cleanIncomplete(boolean dryRun)
}

class CleanupReport {
    int entriesDeleted
    int entriesScanned
    long bytesFreed
    List<String> warnings
    List<HashCode> deletedHashes
}
```

Predicate-driven, ref-count-aware bulk eviction returning a structured `CleanupReport`. Implementations may interact with ref counts (e.g. `DefaultCache.cleanAll` decrements via `CacheDB.removeTaskEntry` and only deletes workdirs when ref hits 0; `GlobalCache.cleanAll` drives the two-phase service protocol).

#### `CacheDropper` — unconditional destructive ops

```groovy
interface CacheDropper {
    /** Drop the per-run index file (implements `nextflow log -clear <runName>`). */
    void deleteIndex()

    /** Drop the entire cache for this uniqueId (implements `nextflow clean -f` and the
      * post-condition of `clean -but` / `clean -after`). */
    void drop()
}
```

Blanket destructive operations — no predicate, no `CleanupReport`, no ref-count check. Separate from `CacheCleaner` because a plugin may reasonably support one without the other.

#### `CacheCleanup` — composite for `openForClean()` return type

```groovy
interface CacheCleanup extends CacheCleaner, CacheDropper, Closeable {}
```

A trivial composite marker. `Cache.openForClean()` returns `CacheCleanup` so the CLI gets a single typed handle covering both selective cleanup and destructive operations without explicit `instanceof` at the entry point. A plugin opts into the composite when it can support both halves; if it supports only one, it implements only `CacheCleaner` or only `CacheDropper` and `openForClean()` throws `UnsupportedOperationException`.

#### Dispatch shape (CLI)

| Verb | Entry point | Required capability |
|------|-------------|---------------------|
| `nextflow log <runName>` | `openForRead()` → `CacheReader` | `CacheReader` |
| `nextflow log -clear <runName>` | `openForClean()` → `CacheCleanup` | `CacheDropper` (via `deleteIndex()`) |
| `nextflow clean -older-than <age>` / `-but <X>` / `-after <X>` | `openForClean()` → `CacheCleanup` | `CacheCleaner` |
| `nextflow clean -f` | `openForClean()` → `CacheCleanup` | `CacheDropper` (via `drop()`) |

When a required capability is missing, `openForRead()` / `openForClean()` throws and the CLI prints `"<impl> does not support <verb>"` (non-zero exit). Silent no-ops would be worse than refusing.

The default implementation (LevelDB / `nf-cloudcache`) implements every capability. Today's `nextflow log` and `nextflow clean` semantics carry over essentially unchanged.

### Section 7 — Plugin API Surface

Plugins extend two `META-INF/services` entry points:

```
META-INF/services/nextflow.cache.CacheFactory
META-INF/services/nextflow.processor.hash.TaskHasherFactory
```

```groovy
abstract class CacheFactory implements ExtensionPoint {

    /** UNCHANGED from pre-pluggable-cache. Every existing CacheFactory subclass
      * (DefaultCacheFactory, CloudCacheFactory, etc.) continues to override only this
      * method. */
    protected abstract CacheDB newInstance(UUID uniqueId, String runName, Path home=null)

    /** NEW — default body wraps the legacy CacheDB in a DefaultCache. Plugins that
      * implement a fully bespoke pluggable Cache override this method directly. */
    Cache create(UUID uniqueId, String runName, Path basePath) {
        return new DefaultCache(newInstance(uniqueId, runName, basePath))
    }

    /** Name of the TaskHasher this cache prefers, resolved via TaskHasherFactory registry */
    String getDefaultTaskHasher() { return 'v2' }

    /** The cache scheme/type for diagnostics (e.g. "default", "cloud", "global"). */
    String getName() { return 'default' }

    /** UNCHANGED static dispatch returning the legacy CacheDB. Used by remaining
      * legacy call sites until they migrate to {@link #createCache}. */
    static CacheDB create(UUID uniqueId, String runName, Path home=null) {
        final all = Plugins.getPriorityExtensions(CacheFactory)
        if( !all ) throw new IllegalStateException("Unable to find Nextflow cache factory")
        return all.first().newInstance(uniqueId, runName, home)
    }

    /** NEW — static dispatch returning the pluggable Cache. */
    static Cache createCache(UUID uniqueId, String runName, Path basePath=null) {
        final all = Plugins.getPriorityExtensions(CacheFactory)
        if( !all ) throw new IllegalStateException("Unable to find Nextflow cache factory")
        return all.first().create(uniqueId, runName, basePath)
    }
}

abstract class TaskHasherFactory implements ExtensionPoint {
    abstract String getName()
    abstract TaskHasher create(TaskRun task)
}
```

**Plugin authoring guidance:**

- A plugin built on the legacy `CacheStore` / `CacheDB` primitives (e.g. `nf-cloudcache`) implements only `newInstance(...)`. The default `create(...)` wraps the returned `CacheDB` in `DefaultCache`; the plugin automatically gets `Cache`, `CacheReader`, and `CacheCleanup` via that wrapping. No additional changes needed.
- A plugin with a fully bespoke Cache implementation (e.g. `nf-global-cache`) overrides `create(...)` and returns its own `Cache` instance. It must still provide `newInstance(...)` (the method is abstract); the typical pattern is to throw `UnsupportedOperationException("This factory does not provide a CacheDB; use create(...) instead.")`.

**Selection at session init**:

```
1. Scan classpath for all CacheFactory implementations.
2. Match config `cache.type` (or default) → select one CacheFactory.
3. Read selected factory's getDefaultTaskHasher() → name.
4. If session-level taskHasher config set, override.
5. If env var NXF_TASK_HASH_VER set, override.
6. Look up named TaskHasherFactory; cache resolved factory on Session.
```

`Session.cache` is typed `Cache` (no second `CacheDB` field). The existing `CacheDB` class remains a public API used internally by `DefaultCache`; alternative `Cache` implementations need not depend on it.

#### `CacheDB` declares the read + drop capabilities

`CacheDB` (public class, signature unchanged) adds `implements Closeable, CacheReader, CacheDropper`. The corresponding method bodies (`getTraceRecord`, `findTraceRecord`, `eachRecord`, `deleteIndex`, `drop`, `close`) already exist verbatim on `CacheDB`; the change is declarative.

`CacheDB.removeTaskEntry` (ref-counted decrement-and-maybe-delete) stays public for backward compatibility but is **not** promoted to any new interface. It is an internal primitive of `DefaultCache.cleanAll`/`cleanByHash`/`cleanOlderThan`/`cleanIncomplete`, together with the dryRun simulation logic that today lives in `CmdClean.wouldRemove` (which moves into `DefaultCache`).

#### Caller migrations

| Caller | Today | After |
|---|---|---|
| `Session.cache` | dual fields: `CacheDB cache` + `Cache cacheImpl` | single `Cache cache`; `Session.start()` uses `CacheFactory.createCache(...).open()` |
| `Session.cleanupAfterRunOK` | manual `eachRecord` + `removeTaskEntry` + workdir delete loop against `CacheDB` | `try (CacheCleanup c = cache.openForClean()) { c.cleanAll(false) }` |
| `Session.notifyTaskSubmit/Complete/Cached` | direct `cache.putIndexAsync` / `putTaskAsync` / `cacheTaskAsync` calls | no cache call — writes owned by `DefaultCache.beginTask` (CLAIM_GRANTED / CACHE_HIT branches) and `DefaultCache.putTaskCacheEntry` |
| `CacheBase.cacheFor(entry)` returning `CacheDB` | single method | split into `readerFor(entry): CacheReader` and `cleanerFor(entry): CacheCleanup` |
| `CmdLog` | `cacheFor(entry).eachRecord(...)` | `try (CacheReader r = readerFor(entry)) { r.eachRecord(...) }`; `-clear` flag uses `cleanerFor(entry).deleteIndex()` |
| `CmdClean` | manual loop calling `cache.removeTaskEntry` per entry + `deleteIndex` + `drop` + the `wouldRemove`/`dryHash` dryRun simulation | `try (CacheCleanup c = cleanerFor(entry)) { c.cleanAll(dryRun); if (!dryRun) { c.deleteIndex(); ...; c.drop() } }`. The dryRun simulation moves into `DefaultCache.cleanAll`. |
| `WaveDebugCmd` | `CacheFactory.create(...)` returning `CacheDB` | `CacheFactory.createCache(...).openForRead()` returning `CacheReader`; same method calls (`findTraceRecord`/`getTraceRecord`), narrower receiver type |

### Section 8 — Backward Compatibility & Migration

With no plugin installed and no opt-in config, every behavior is byte-identical to today.

**Default behavior preservation per seam:**

| Seam | Default impl behavior | Risk vs today |
|------|----------------------|---------------|
| `TaskHasher` | `v2` selected by factory; `compute()` produces same hash for same task | None — interface extracted from concrete class without changing logic |
| `beginTask` | Today's loop relocated: load `TaskEntry`, `LockManager.acquire` + workdir.mkdirs | Low — pure relocation |
| `endTask` | No-op for ref counting; lock already released | None — additive |
| `getTaskCacheEntry` | Loads `TaskEntry`, scans workdir for `.exitcode`/`.command.out`/output files; assembles `OutputBinding` map | Low — pure relocation. The scan moves from `TaskProcessor` into the cache impl; total filesystem work is identical |
| `putTaskCacheEntry` | Same as today's `writeTaskEntry0`/`writeTaskIndex0` (TaskEntry serialized; output bindings discarded) | None — same on-disk format |
| `adopt` | Returns input bindings unchanged + `KEEP` | None — additive no-op |
| `notifyPublish` / `notifyFilePort` | Empty methods | None — additive no-op |

**On-disk format guarantees:**

- LevelDB cache files: unchanged.
- `nf-cloudcache` file layout: unchanged (`uniqueId/<hash>` files + `index.<runName>`).
- Existing caches from prior Nextflow versions: continue to be readable.

**Code paths NOT in scope for behavior change:**

- `checkStoredOutput` / `storeDir` semantics — bypasses the cache entirely; no cache calls.
- Error handling / retry strategies (`maxErrors`, `maxRetries`, dynamic resources) — `beginTask`/`endTask` wrap, do not replace, the existing logic.
- `TraceObserverV2` events — independent.
- Lineage tracking (`nf-lineage`) — independent.

**Opt-in model:**

```groovy
// nextflow.config — current behavior, no opt-in needed (defaults preserve today)

// Opt in to alternative cache implementation (e.g., global cache plugin)
plugins {
    id 'nf-global-cache@1.0.0'   // hypothetical plugin name
}
cache {
    type = 'global'              // selects the plugin's CacheFactory
}
```

The hasher is selected automatically based on `CacheFactory.getDefaultTaskHasher()` declared by the active factory.

**Phased rollout:**

1. **Phase A — foundation**: ship `TaskHasher` interface (PR #6927) + plugin-extensible factory. Zero behavior change.
2. **Phase B — cache contract**: introduce `beginTask`/`endTask`/`getTaskCacheEntry`/`putTaskCacheEntry` with default implementations wrapping today's logic. Existing `getTaskEntry` method kept as a transitional alias (deprecated). Zero behavior change.
3. **Phase C — adoption + events**: add `adopt`, `notifyPublish`, `notifyFilePort`. Zero behavior change.
4. **Phase D — capability interfaces**: add `CacheReader`, `CacheCleaner`, `CacheDropper`, `CacheCleanup`. Route `nextflow log` through `openForRead()` / `CacheReader` and `nextflow cache clean` through `openForClean()` / `CacheCleanup`. Migrate `CacheBase`, `CmdLog`, `CmdClean`, `Session.cleanupAfterRunOK`, `WaveDebugCmd`. Zero behavior change.
5. **Phase E — validation**: a real alternative cache implementation (the global cache plugin) validates the contract end-to-end before declaring it stable.

Each phase ships independently; no big-bang merge.

### Section 9 — Risks

1. **Event pairing bugs (`beginTask`/`endTask` not paired in some error path)** — high likelihood × high impact. A single missed `endTask` permanently leaks a ref count, pinning a file forever in a ref-counted cache. Mitigation: dedicated test class with parameterized failure injection across every error strategy (`retry`, `terminate`, `ignore`, `finish`); audit every `return`/`throw` in `TaskProcessor` between the call sites; consider a try-with-resources style helper that guarantees `endTask` fires.
2. **Adoption-vs-DELETE invariant violations** — medium likelihood × high impact. If `adopt()` returns `DELETE` but does not rewrite all FILE bindings, downstream consumers see vanished files. Mitigation: enforce as a runtime contract — assertion fails the run loudly at adopt time; comprehensive test for the global-cache adoption path.
3. **Hash bump policy drift between cache implementations** — low likelihood × medium impact. Default uses `bump(hash, tries) = hash + tries`. If a global-cache implementation uses a different bump function, two pipelines hitting the same task with different cache configs produce different hash sequences. Mitigation: document the bump function as a contract concern in plugin authoring docs.

**`collectOutputs` semantics drift** is structurally prevented: visibility on `TaskProcessor.collectOutputs(TaskRun)` is bumped from `protected` to `public`, and the default cache implementation calls `task.processor.collectOutputs(task)` directly. Both the success path and the cache-restore path exercise the same single method.

### Known Limitations

- **Single hasher per session.** Mixed pipelines cannot have different hashers per process. Future extension if needed.
- **Adoption is sync, not background.** A slow `adopt()` (large output upload) blocks task finalization. Async adoption complicates failure semantics ("upload failed after channel binding") and is deferred until a real bottleneck is observed.
- **`storeDir` outputs bypass the cache.** Consistent with today.
- **Workflow-level outputs only emit `notifyPublish`.** Cache does not ref-count workflow output state — that is the workflow's permanent product, outside cache lifecycle.
- **Failed tasks do not adopt.** Partial outputs from a failed run never enter the cache.
- **Trace observers are independent.** Cache implementations that also want trace data still consume `TraceObserverV2` events separately.
- **Cleanup-during-execution remains a cache-impl concern.** With authoritative-via-events the implementation has enough information not to delete in-use files; correctness is on the implementation, not the contract.

### Out of Scope (deferred to follow-up specs)

- Per-process hasher selection.
- Async / batched event emission.
- Cross-pipeline lineage integration.
- Performance tuning specifics (bloom filters, fast-hash → deep-hash cache from `adr/20260202-global-cache.md` Alternative 7).
- Specific global cache implementation details — covered by `adr/20260202-global-cache.md` and any plugin-specific specs.

## Acceptance Criteria

- All five seams are introduced with default implementations that pass the existing test suite without modification.
- A reference alternative implementation (the global cache plugin) wires all interfaces end-to-end and demonstrates outputs-shaped restore, workdir adoption with `DELETE`, and ref-counted file usage.
- `nextflow cache clean ...` continues to work for the default implementation; refuses gracefully for caches that do not implement the corresponding capability (`CacheCleaner` for selective eviction, `CacheDropper` for `-f` / drop).
- `LockManager` is no longer referenced from `TaskProcessor`.
- `collectOutputs` is no longer called from the cache-restore path in `TaskProcessor`.
- LevelDB and `nf-cloudcache` on-disk formats are unchanged; pre-existing caches resume successfully.
- A failure-injection test demonstrates `beginTask`/`endTask` pairing across every error strategy.

---

## Clarifications

### Session 2026-05-07

- Q: How should `TaskProcessor` handle exceptions thrown by cache plugin methods (`beginTask`, `endTask`, `adopt`, `putTaskCacheEntry`, `notifyPublish`, `notifyFilePort`)? → A: Per-method matrix anchored to today's behavior — `beginTask` exception ⇒ log warning, treat as cache miss, proceed to claim+execute (preserves today's `getTaskEntry`/`checkCachedOutput` catch-and-fall-through at TaskProcessor.groovy:810–823); `endTask` ⇒ log warning, continue (cache plugin must tolerate missed `endTask` via its own TTL/recovery for crash cases anyway); `adopt` ⇒ fail the task and honor `errorStrategy` (cache integration failed, downstream URIs may be invalid, safer to fail loudly); `putTaskCacheEntry` ⇒ log warning, continue (matches today's async-swallow behavior; the task already succeeded); `notifyPublish` / `notifyFilePort` ⇒ log warning, continue (best-effort observability).
- Q: How do existing `CacheDB` consumers and `CacheStore` plugins migrate to the expanded interface? → A: Introduce a new public interface `nextflow.cache.Cache` carrying `beginTask`, `endTask`, `adopt`, `getTaskCacheEntry`, `putTaskCacheEntry`, `notifyPublish`, `notifyFilePort`. `TaskProcessor` is retyped to use `Cache` as a single code path (no `instanceof` checks). The existing `CacheDB` class is left untouched (no rename, no `@Deprecated` annotation) and continues as an internal building block of the default `Cache` implementation, which wraps a `CacheDB` (which itself wraps a `CacheStore`). Existing `CacheStore`-based plugins (e.g., `CloudCacheStore` in `nf-cloudcache`) keep working unchanged — they are picked up by the default `Cache` implementation's wrapper. Plugins that need full extensibility implement `Cache` directly, bypassing both `CacheDB` and `CacheStore`. Methods that exist only on `CacheDB` today (`eachRecord`, `findTraceRecord`, `incTaskEntry`, `removeTaskEntry`, `cacheTaskAsync`) stay there for the default impl's internal use and any direct access (CLI commands, observers).
- Q: How are `TaskHasher` factory name conflicts resolved? → A: `v1` and `v2` are reserved built-in names; plugins MUST NOT register a `TaskHasherFactory` with either name. Attempting to do so aborts session start with a clear error identifying the offending plugin. Plugin-vs-plugin name collisions (two non-built-in factories registering the same name) likewise abort session start with an error listing both plugins. Hash semantics are load-bearing for cache correctness — silently picking one of two "global" implementations could produce wrong hashes and cause silent cache misses or false hits if the two implementations happen to share a hash space. Plugin authors are expected to namespace their factory names (e.g. `nf-global-cache:global` rather than just `global`).
- Q: What observability does Nextflow core provide around `Cache` dispatch? → A: Preserve today's user-facing log lines unchanged (`info` cache-hit announcement at `TaskProcessor.groovy:991`, `trace` lookup details at `TaskProcessor.groovy:816`). Add a single `debug`-level structured log line per `Cache` method dispatch — for example: `[abc1234] Cache.beginTask outcome=CACHE_HIT finalHash=def5678` and `[abc1234] Cache.adopt disposition=DELETE bindings_rewritten=3`. Cache-plugin failures (per the exception-handling clarification above) continue to be logged at `warn`. No core-side metrics, counters, or `TraceObserverV2`-style cache events are introduced; cache plugins are responsible for their own metrics, dashboards, and ops integrations. Operators get enough diagnostic surface in the Nextflow log to file actionable bug reports against plugin authors; plugin authors get freedom to ship richer telemetry without core changes.
- Q: What lifetime guarantees does `TaskCacheEntry.stdoutUri` carry? → A: `stdoutUri` MUST be readable for the lifetime of the cache entry. Cache plugins MUST replicate stdout into managed storage (i.e. outside the workdir) before returning `WorkdirDisposition.DELETE` from `adopt`, the same way they MUST rewrite all FILE-typed bindings. The default implementation returning `KEEP` retains stdout in the workdir as today. This makes stdout a load-bearing artifact for trace records, `--with-trace`, and operator-facing error reports — operators can rely on cached stdout being available regardless of which cache implementation is active.

---

## User Scenarios & Testing

### User Story 1 — Existing pipeline runs unchanged (Priority: P1)

A pipeline operator who has not installed any cache plugin and not changed any cache config runs a pipeline with `-resume`. Cache hits, cache misses, retries, and failures behave exactly as they did before the refactor.

**Why this priority**: Backward compatibility is the load-bearing constraint of this entire refactor. A regression here breaks every existing Nextflow user. Without this, no other story is allowed to ship.

**Independent Test**: Run the full existing Nextflow test suite (unit + integration + cloud validation) against a build with the refactor; verify zero failures and zero behavioral diffs versus master. Run a real pipeline (e.g., `nf-core/sarek`) with `-resume` from a pre-refactor cache; verify all tasks resume as cached.

**Acceptance Scenarios**:

1. **Given** a pipeline previously executed without the refactor (LevelDB cache on disk), **When** the operator runs the same pipeline on a refactored Nextflow with `-resume`, **Then** every task that would have been cached is restored from cache with the same hash and the same output channels.
2. **Given** a pipeline using `nf-cloudcache` on S3, **When** the operator runs `-resume` after the refactor, **Then** cached entries written before the refactor are restored without modification.
3. **Given** a task that fails and is configured with `errorStrategy 'retry'`, **When** the task fails on attempt 1 and succeeds on attempt 2, **Then** the retry produces the same hash sequence as today and stores its result in cache identically.
4. **Given** a task whose hash is computed from inputs, container, and script, **When** the operator inspects the hash on a refactored build, **Then** the hash bytes are identical to those produced on master for the same task.

---

### User Story 2 — Plugin author ships an alternative cache implementation (Priority: P1)

A plugin author creates a Nextflow plugin that provides a custom `Cache` and `TaskHasher` (for example, a global cache backed by a service or a different storage strategy). The plugin can override hash computation, cache lookup, claim coordination, output adoption, and file-usage tracking without modifying core Nextflow code.

**Why this priority**: This is the entire point of the refactor. Without it the global cache and any future cache variant cannot ship.

**Independent Test**: Build a minimal sample plugin that registers a `CacheFactory` whose `Cache` implementation logs every method call; run a small pipeline with the plugin enabled; verify all expected lifecycle methods (`beginTask`, `adopt`, `putTaskCacheEntry`, `endTask`, `notifyPublish`, `notifyFilePort`) fire in the correct order with correct arguments.

**Acceptance Scenarios**:

1. **Given** a plugin registers a `CacheFactory` with name `"custom"`, **When** the operator sets `cache.type = 'custom'` in `nextflow.config`, **Then** Nextflow resolves the active cache via the plugin's factory and routes all cache-related calls to its `Cache` implementation.
2. **Given** the active `CacheFactory` declares `getDefaultTaskHasher() = "custom"` and a `TaskHasherFactory` named `"custom"` is registered, **When** the session initialises, **Then** `Session.hashStrategy` resolves to the custom hasher without requiring user config.
3. **Given** a custom `Cache.beginTask` returns `CACHE_HIT` with a `TaskCacheEntry`, **When** `TaskProcessor` consumes the result, **Then** outputs bind to channels using only the entry's contents, with no workdir scan performed in core.
4. **Given** a custom `Cache.adopt` returns `WorkdirDisposition.DELETE` and rewrites all FILE bindings to URIs outside the workdir, **When** `finalizeTask` completes, **Then** the workdir is asynchronously deleted and downstream tasks consume the rewritten URIs.
5. **Given** a plugin author wants to ship the implementation, **When** they package it as a standard Nextflow plugin with `META-INF/services/nextflow.cache.CacheFactory` and `META-INF/services/nextflow.processor.hash.TaskHasherFactory` entries, **Then** the plugin loads via the existing plugin mechanism with no core code changes required.

---

### User Story 3 — Cache plugin tracks file usage for safe cleanup (Priority: P2)

A plugin author implementing reference-counted cleanup uses the `notifyPublish` and `notifyFilePort` events plus the inputs argument of `beginTask` and the success flag of `endTask` to maintain authoritative ref counts. The plugin can determine when a file is in use and never delete a referenced file during execution.

**Why this priority**: Cleanup correctness is the second motivating use case after sharing. Authoritative-via-events is the contract we commit to so that plugins can build correct cleanup without polling Nextflow.

**Independent Test**: Run a pipeline with a sample plugin that records every `beginTask`/`endTask`/`notifyPublish`/`notifyFilePort` to an event log. After the run, verify: every `beginTask` is paired with exactly one `endTask`; no `endTask` precedes the channel-binding step that would create a downstream consumer; `notifyPublish` fires for every published file.

**Acceptance Scenarios**:

1. **Given** task A produces output O consumed by task B, **When** the run executes with both tasks, **Then** the cache plugin observes `beginTask(A)` → `adopt(A, …)` → `putTaskCacheEntry(A, …)` → channel binding → `endTask(A, success=true)` → `beginTask(B, inputs=[O, …])` → … → `endTask(B, …)`, in that order.
2. **Given** a task fails with `errorStrategy 'terminate'`, **When** the failure propagates, **Then** the cache plugin observes `endTask(task, success=false)` exactly once for the failed task.
3. **Given** a `publishDir` directive, **When** an output file is successfully published, **Then** the cache plugin observes `notifyPublish(source, destination, mode)` exactly once with the actual mode used (COPY, SYMLINK, MOVE, or LINK).
4. **Given** an input file at a remote URI, **When** `FilePorter` downloads it locally for staging, **Then** the cache plugin observes `notifyFilePort(remote, local, task)` after the download completes.

---

### User Story 4 — Operator opts into an alternative cache implementation (Priority: P2)

A pipeline operator wants to use an alternative cache (for example, the global cache plugin). They install the plugin, set the cache type in config, and run their pipeline. Their cached outputs land in the cache plugin's managed storage; the workdir becomes disposable; reruns of the same task across different sessions hit the alternative cache.

**Why this priority**: The user-facing story for non-default caching. Important but predicated on plugin authors having shipped (Story 2).

**Independent Test**: With a sample alternative cache plugin installed, run a small pipeline; verify outputs land in plugin-managed storage; rerun the same pipeline from a clean workdir base; verify task hashes hit the cache and no tasks re-execute.

**Acceptance Scenarios**:

1. **Given** the operator adds the plugin to `nextflow.config` and sets `cache.type` to the plugin's name, **When** they run a pipeline, **Then** Nextflow loads the plugin, dispatches all cache calls to it, and uses its declared default hasher.
2. **Given** the plugin's `adopt` returns `DELETE`, **When** a task succeeds, **Then** the workdir is deleted after finalization and downstream tasks read from the plugin's URIs.
3. **Given** the operator removes the plugin and reverts the config, **When** they run the pipeline again, **Then** Nextflow reverts to the default cache implementation without errors.

---

### User Story 5 — Operator runs `nextflow cache clean` (Priority: P3)

A pipeline operator runs `nextflow cache clean` (or a variant: `-older-than`, `-hash`, `-incomplete`, `-dry-run`) against a pipeline. The CLI calls `Cache.openForClean()` to obtain a `CacheCleanup` handle and dispatches to the right method. The default cache cleans as today; an alternative cache implementation that does not support cleanup throws `UnsupportedOperationException` from `openForClean()` and the CLI refuses gracefully with a clear message.

**Why this priority**: Operational story for users with cleanup needs. Lower priority than the architectural seams because existing `nextflow clean` semantics already work for the default case; the new behavior is the refusal path for caches that do not implement `CacheCleaner` / `CacheDropper`.

**Independent Test**: Run `nextflow cache clean -dry-run` against the default cache; verify entries are reported. Repeat against a cache plugin whose `openForClean()` throws; verify the CLI exits non-zero with a clear message and no entries are touched.

**Acceptance Scenarios**:

1. **Given** the default cache, **When** the operator runs `nextflow cache clean -older-than 7d`, **Then** entries older than seven days are deleted and a summary report is printed.
2. **Given** an alternative cache plugin whose `openForClean()` returns a `CacheCleanup`, **When** the operator runs `nextflow cache clean`, **Then** the CLI dispatches to the plugin's implementation and reports its results.
3. **Given** an alternative cache plugin whose `openForClean()` throws `UnsupportedOperationException`, **When** the operator runs `nextflow cache clean`, **Then** the CLI prints "this cache implementation does not support manual cleanup" and exits with non-zero status.
4. **Given** any cache implementation, **When** the operator runs `nextflow cache clean -dry-run`, **Then** the CLI reports what would be deleted without modifying any state.

---

### Edge Cases

- **Task fails between `beginTask` and `endTask`**: `endTask(success=false)` MUST still fire from the appropriate finalization path. Audit covers all `errorStrategy` values (retry, terminate, ignore, finish) and dynamic-resources retry.
- **JVM crashes mid-task**: Documented gap. Cache implementation may include crash-recovery (e.g., stale `.lock` markers cleaned by age) but core makes no on-disk transactionality guarantee.
- **Cache plugin returns `WorkdirDisposition.DELETE` but leaves a FILE binding pointing into the workdir**: Runtime assertion fails the run loudly at adopt time.
- **Plugin registers a `TaskHasherFactory` with a reserved built-in name (`v1` or `v2`)**: Session start aborts with a clear error identifying the offending plugin.
- **Two plugins register `TaskHasherFactory` with the same non-built-in name**: Session start aborts with an error listing both conflicting plugins. Plugin authors are expected to namespace their factory names (e.g. `nf-global-cache:global`).
- **Workdir delete fails after `disposition=DELETE`**: Logged but does not fail the task or the workflow; the cache entry is already persisted with rewritten URIs and downstream tasks are unaffected.
- **`beginTask` is called with `tryHit=false` (retry path)**: Cache implementation skips hit lookup and proceeds directly to claim acquisition with bumped hash.
- **Concurrent `beginTask` calls within the same session for the same hash**: Default implementation serialises via `LockManager`; alternative implementations are responsible for their own intra-session and cross-session serialisation.
- **Cache plugin hangs in `beginTask` or `adopt`**: No core-side timeout. Plugin authors are responsible for their own timeouts and circuit breakers.
- **Cache plugin throws an exception from `beginTask`**: Log warning, treat as cache miss, fall through to claim+execute (preserves today's behavior at `TaskProcessor.groovy:810–823`).
- **Cache plugin throws from `endTask`, `putTaskCacheEntry`, `notifyPublish`, or `notifyFilePort`**: Log warning and continue. These events are best-effort; the task itself is unaffected.
- **Cache plugin throws from `adopt`**: Fail the task and honor the configured `errorStrategy` (`retry`, `terminate`, `ignore`, `finish`). Outputs were produced but cache integration failed and downstream URIs may be invalid; failing loudly is safer than masking.
- **`storeDir`-based outputs**: Bypass the cache entirely — `checkStoredOutput` does not call `beginTask`/`endTask`/`adopt`/`putTaskCacheEntry`.
- **Failed task adoption**: `adopt` is never called for failed tasks; partial outputs from a failed run never enter the cache.

## Requirements

### Functional Requirements

- **FR-001**: System MUST select the active `TaskHasher` via a plugin-extensible factory at session initialisation, resolving by name in the precedence order: cache-driven default, session config, env var, hardcoded `v2`.
- **FR-002**: With no plugin registered and no cache config set, the resolved `TaskHasher` MUST produce hash bytes byte-identical to those produced by master for the same task.
- **FR-003**: System MUST call `cache.beginTask(BeginTaskRequest)` exactly once before each task's execution decision (cache hit or claim acquisition).
- **FR-004**: System MUST call `cache.endTask(EndTaskRequest)` exactly once per `beginTask`, in both success and failure paths, after channel binding completes.
- **FR-005**: System MUST call `cache.adopt(task, outputs)` exactly once for each successfully-completed task, before the cache entry is persisted and before output channels are bound.
- **FR-006**: System MUST honour `WorkdirDisposition.DELETE` by deleting the task workdir asynchronously after task finalisation.
- **FR-007**: System MUST validate that when `adopt` returns `DELETE`, every FILE-typed `OutputBinding` and the resulting `TaskCacheEntry.stdoutUri` resolve to URIs outside the workdir; violations MUST fail the run.
- **FR-008**: System MUST emit `notifyPublish(source, destination, mode)` exactly once after every successful publish operation performed by `PublishDir` or `PublishOp`.
- **FR-009**: System MUST emit `notifyFilePort(remote, local, task)` exactly once after every successful download performed by `FilePorter`.
- **FR-010**: `cache.getTaskCacheEntry(hash, processor)` MUST return enough information to bind output channels and finalize the task without any workdir scan in `TaskProcessor`.
- **FR-011**: System MUST NOT call `collectOutputs` or read `.exitcode`/`.command.out` from `TaskProcessor` on the cache-restore path.
- **FR-012**: `TaskProcessor.collectOutputs(TaskRun)` MUST be the single source of truth for output collection semantics; both the success path and the cache-restore path (in the default cache implementation) MUST invoke it directly without duplicating logic.
- **FR-013**: Default cache implementation MUST preserve existing on-disk format for both LevelDB and `nf-cloudcache`; cache files written by prior Nextflow versions MUST be readable.
- **FR-014**: System MUST dispatch `nextflow cache clean` to the active cache when `Cache.openForClean()` returns a `CacheCleanup` handle.
- **FR-015**: System MUST refuse `nextflow cache clean` with a clear message and non-zero exit when the active cache's `openForClean()` throws `UnsupportedOperationException`.
- **FR-016**: `LockManager` MUST NOT be referenced from `TaskProcessor` after the refactor; intra-session serialisation MUST be an internal concern of the default cache implementation.
- **FR-017**: Cache implementations MUST own the hash-bump policy used for collision avoidance; `TaskProcessor` MUST NOT contain a retry/bump loop for cache-related collisions.
- **FR-018**: Plugin authors MUST be able to register a `CacheFactory` and a `TaskHasherFactory` via standard `META-INF/services` discovery without modifying core code.
- **FR-019**: When `beginTask` returns `CACHE_HIT`, the cache implementation is responsible for ref-counting any reused outputs internally; no separate "cache hit" notification is emitted to the cache.
- **FR-020**: System MUST NOT call `adopt` or `putTaskCacheEntry` for tasks that fail or for tasks served via `storeDir`.
- **FR-021**: When a cache plugin method throws an exception, the system MUST behave as follows: `beginTask` ⇒ log warning and proceed as cache miss; `endTask`, `putTaskCacheEntry`, `notifyPublish`, `notifyFilePort` ⇒ log warning and continue; `adopt` ⇒ fail the task and honor the configured `errorStrategy`. The cache plugin's exception MUST NOT abort the workflow except as a normal task-failure consequence of `errorStrategy 'terminate'`.
- **FR-022**: `v1` and `v2` are reserved `TaskHasher` factory names. The system MUST abort session start with a clear error if a plugin registers a `TaskHasherFactory` whose `getName()` returns a reserved name.
- **FR-023**: The system MUST abort session start with a clear error if two `TaskHasherFactory` registrations resolve to the same non-built-in name; the error MUST identify both conflicting plugins.
- **FR-024**: The system MUST emit one `debug`-level structured log line per `Cache` method dispatch including the task hash, method name, and a method-specific summary (outcome for `beginTask`; disposition + bindings-rewritten count for `adopt`; success flag for `endTask`; etc.). User-facing logs at `info` and `trace` MUST remain unchanged from today's output.
- **FR-025**: Core MUST NOT introduce new metrics, counters, or `TraceObserverV2` events for `Cache` dispatch. Cache plugins are responsible for their own metrics and ops integrations.
- **FR-026**: `TaskCacheEntry.stdoutUri` MUST be readable for the lifetime of the cache entry. Cache implementations that return `WorkdirDisposition.DELETE` from `adopt` MUST replicate stdout into managed storage outside the workdir before returning.

### Key Entities

- **TaskHasher**: Computes the content-derived hash for a task. Pluggable per session.
- **TaskHasherFactory**: Plugin extension point that produces a named `TaskHasher` implementation.
- **Cache**: New public interface that owns hash-keyed lookup, claim coordination, output adoption, persistence, and file-usage events. Pluggable per session via `CacheFactory`.
- **CacheFactory**: Plugin extension point that produces a named `Cache` and declares a preferred `TaskHasher` name.
- **CacheDB**: Existing concrete class (untouched) that wraps a `CacheStore` with async writer and Kryo serialisation. Used internally as a building block of `DefaultCache`; not part of the public plugin contract going forward.
- **DefaultCache**: Default `Cache` implementation. Wraps a `CacheDB` (which wraps a `CacheStore`) plus a `LockManager`. Preserves today's behavior byte-for-byte and on-disk format unchanged.
- **TaskCacheEntry**: Outputs-shaped restore artifact returned by the cache. Carries trace, context, exit code, stdout URI, and an `OutputBinding` map sufficient to bind channels without a workdir scan.
- **OutputBinding**: Per-output-parameter binding (FILE / VALUE / EVAL) with value(s) and optional metadata.
- **CacheResolution**: Outcome of `beginTask` — `CACHE_HIT` (with entry) or `CLAIM_GRANTED` (with finalHash).
- **AdoptionResult**: Outcome of `adopt` — possibly-rewritten output bindings plus `WorkdirDisposition` (KEEP/DELETE).
- **CacheReader**: Optional capability interface (`getTraceRecord` / `findTraceRecord` / `eachRecord`) returned by `Cache.openForRead()`. Enables `nextflow log` and per-run reporting.
- **CacheCleaner**: Optional capability interface for predicate-driven, ref-count-aware bulk eviction (`cleanAll` / `cleanByHash` / `cleanOlderThan` / `cleanIncomplete`).
- **CacheDropper**: Optional capability interface for unconditional destructive operations (`deleteIndex` / `drop`).
- **CacheCleanup**: Composite of `CacheCleaner` + `CacheDropper` + `Closeable`; the narrowed return type of `Cache.openForClean()`. Enables `nextflow cache clean` dispatch.

## Success Criteria

### Measurable Outcomes

- **SC-001**: Pre-refactor caches resume with 100% of cached tasks restored on a refactored build (zero re-execution induced by the refactor itself).
- **SC-002**: The full existing Nextflow test suite (unit, integration, cloud validation) passes with zero new failures on a refactored build with no plugins enabled.
- **SC-003**: Hash bytes produced by the default `TaskHasher` (`v2`) are byte-identical to master for 100% of tasks across a representative pipeline corpus (e.g., `nf-core/sarek`, `nf-core/rnaseq`).
- **SC-004**: A reference alternative cache plugin demonstrates outputs-shaped restore, workdir adoption with `DELETE`, and ref-counted file usage end-to-end on a pipeline of at least 100 tasks with no core code modifications.
- **SC-005**: Failure-injection tests pass for every supported `errorStrategy` (`retry`, `terminate`, `ignore`, `finish`) and for dynamic-resources retry, demonstrating exactly one `endTask` call per `beginTask` call across all paths.
- **SC-006**: Plugin authors can ship a working `Cache` implementation as a standalone Nextflow plugin without changes to core (validated by the reference implementation shipping in a separate repository or plugin).
- **SC-007**: `nextflow cache clean` continues to work for the default implementation with the same UX as today; refuses gracefully (non-zero exit, clear message, zero state changes) for cache plugins whose `openForClean()` throws.
- **SC-008**: After enabling an alternative cache plugin that returns `DELETE`, workdirs are removed at a rate of at least 95% of completed tasks (the remaining 5% accounts for rare async-delete failures, which MUST be logged and MUST NOT block the workflow).

## Assumptions

- Plugins are discovered via the existing `pf4j` mechanism with `META-INF/services` registration; no new plugin discovery infrastructure is required.
- A single `CacheFactory` is selected per session; mixed cache configurations across processes are out of scope.
- A single `TaskHasher` is selected per session (per-process selection deferred to a follow-up spec).
- The reference alternative cache plugin used for end-to-end validation is the global cache plugin tracked in `adr/20260202-global-cache.md`; the spec does not depend on any specific cache implementation other than the default.
- `TraceObserverV2` and `nf-lineage` remain independent of the cache contract; cache plugins that also need trace or lineage data subscribe to them separately.
- `storeDir`-based output handling is unchanged; tasks served via `storeDir` continue to bypass the cache entirely.
