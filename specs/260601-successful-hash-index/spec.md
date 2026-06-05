# Successful-Hash Index for the Default Cache

- Author: Jorge Ejarque
- Date: 2026-06-01
- Status: draft

## Summary

Add a small index to Nextflow's default task cache that maps a task's
**content hash** to the **final hash** of its first successful execution. On
resume, the cache resolves a completed task in a single index lookup instead of
walking the sequence of per-attempt hashes. This is an alternative to the
`previousTryCount` approach proposed in PR #6903 for the same underlying
problem (issue #6884): a task that succeeds only after one or more failed or
aborted attempts must be re-discoverable on resume without re-executing.

## Problem

### How the default cache resolves a task today

When a process invocation is evaluated, Nextflow computes a content-based hash
of the task (`new TaskHasher(task).compute()` in
`TaskProcessor.checkCachedOrLaunchTask`). It then enters an *iterate-and-bump*
loop:

```
int tries = task.failCount + 1
while( true ) {
    hash = H(contentHash + tries)              // per-attempt "final hash"
    entry = session.cache.getTaskEntry(hash)
    resumeDir = entry?.trace?.workDir
    exists = resumeDir?.exists()
    if( shouldTryCache && exists && entry.trace.isCompleted()
            && checkCachedOutput(...) )
        break                                  // cache hit -> resume
    if( exists ) { tries++; continue }         // attempt slot taken -> bump
    // else: allocate workdir for `hash` and launch
}
```

The `tries` counter mixed into the hash gives every retry attempt a distinct
hash and workdir. A first-attempt success lives at `H(contentHash + 1)`; a
success after two failures lives at `H(contentHash + 3)`.

### The failure mode (issue #6884)

A task that fails on attempt 1 and succeeds on attempt 2 has cache entries at
`H(contentHash + 1)` (failed/aborted) and `H(contentHash + 2)` (completed). On
resume, the loop must walk past the failed attempt to find the success. That
walk is only correct as long as every intermediate attempt's workdir and cache
entry still exist on disk and carry the expected status. When that coupling
breaks — a failed-attempt workdir was cleaned, an entry was written without the
expected status, the run was killed mid-retry — the walk can stop early and
**re-execute a task that had already completed successfully**. The reported
reproduction:

1. Configure a task to fail on attempt 1, succeed on attempt 2.
2. Run until that task completes; kill the pipeline.
3. Resume; let it complete; kill again.
4. Resume once more — the already-completed task re-executes instead of
   resuming.

PR #6903 addresses this by recording the count of prior unsuccessful attempts
(`previousTryCount`) on `TaskRun` and folding failed/aborted history into the
starting `tries`, so the walk lands on the right slot. That keeps the fragile
walk but makes its starting point recoverable.

### This proposal

Record the answer directly. When a task succeeds, store
`contentHash -> finalHash` in an index. On resume, one lookup returns the hash
of the successful execution, and the cache validates and resumes that entry
without walking the attempt sequence at all. The walk reverts to its original,
narrow role: allocating a fresh hash slot when launching or retrying.

## Goals

- Resolve a completed task on resume in a single index lookup, regardless of
  how many failed/aborted attempts preceded it.
- Fully replace the `previousTryCount` mechanism for the default cache — one
  mechanism, not two.
- Always on, no configuration flag. A strict improvement with graceful
  fallback to today's behaviour.
- Keep the change small and reviewable, and keep the proven launch/collision
  loop untouched.

## Non-goals

- Cross-pipeline / cross-machine cache sharing, content-based input hashing,
  atomic-create coordination, or ref-counting redesign. Out of scope.
- Migrating existing caches. Old caches simply lack the index and fall back to
  current behaviour; a subsequent run repopulates the index.

## Terminology

- **contentHash** — `new TaskHasher(task).compute()`; the task's
  content-based identity *before* per-attempt try-mixing. The index **key**.
- **finalHash** — `H(contentHash + tries)`; the value assigned to `task.hash`
  and used as the cache-entry key. The index **value**.

## Design

### Approach

A surgical fast-path at the top of `checkCachedOrLaunchTask`. When the index
returns a pointer and the target entry validates, resume directly and return.
Otherwise fall through to the existing iterate-and-bump loop unchanged. The
loop continues to serve two roles: write-side collision avoidance for retries
within a run, and the fallback read path for caches that predate the index. On
any successful resume or execution, the index pointer is written (or rewritten).

To keep `checkCachedOrLaunchTask` readable, the fast-path's lookup/validation
logic is **extracted into a private `resumeFromSuccessfulHashIndex(task)` helper** rather
than inlined; the method itself only gains the content-hash capture and a
single `if( shouldTryCache && resumeFromSuccessfulHashIndex(task) ) return`. The original
loop body is **not** refactored — it stays byte-for-byte as before. This is a
narrower change than splitting the whole method into `tryResumeFromIndex()` +
`launchTask()`, which was rejected as a larger blast radius on a
concurrency-sensitive method with no behavioural benefit.

### Component changes

**1. `CacheStore` interface** (`nextflow.cache.CacheStore`) and both
implementations.

```java
HashCode getSuccessfulHash(HashCode contentHash)            // null if absent
void     putSuccessfulHash(HashCode contentHash, HashCode finalHash)
void     deleteSuccessfulHash(HashCode contentHash)
```

The index is **forward only** (`contentHash → finalHash`): `finalHash =
H(contentHash + tries)` is one-way, and cache entries / the run-index are keyed
by finalHash, so `removeTaskEntry(finalHash)` cannot derive the index key from
the hash alone. To make per-entry cleanup possible, the **content hash is
stored in the cache entry record itself** (a fourth slot — see *CacheDB*
below). When an entry is deleted, `removeTaskEntry` recovers the content hash
from the record and **conditionally** removes the index pointer — only when the
pointer actually targets the entry being deleted (a failed attempt shares the
same content hash, but the pointer aims at the *successful* finalHash, so it
must be preserved). A full cache `drop()` still wipes the index wholesale.

- `DefaultCacheStore` (LevelDB): store the pointer under a **namespaced key** —
  a 1-byte prefix (`0x01`) prepended to `contentHash.asBytes()`, value =
  `finalHash.asBytes()`.

  **Key-space contract.** The single LevelDB instance hosts two disjoint key
  namespaces — task entries at `key = <hash>` (exactly `KEY_SIZE` bytes) and
  index pointers at `key = 0x01 + <hash>` (`KEY_SIZE + 1` bytes). The extra byte
  makes the two lengths different, so an index key can never byte-collide with an
  entry key (LevelDB compares the full byte sequence). This is the idiomatic
  embedded-KV pattern (one ordered store, prefixed keys for multiple logical
  record types — as used by, e.g., Bitcoin Core / go-ethereum, and the role
  RocksDB column families fill); it is preferred over a second LevelDB instance
  here because the index shares the entries' lifecycle (created, cleaned, and
  `drop()`-ed together), access pattern (point lookups), and durability, so none
  of the triggers that justify a separate store apply — and a shared store keeps
  one lifecycle and preserves the option of atomic batched writes.

  *Invariant to preserve:* all key construction goes through a single
  `successfulIndexKey()` builder; the store is **never** iterated by LevelDB key (record
  enumeration uses the separate flat index file), so this heterogeneous keyspace
  is invisible today. Should a LevelDB key scan ever be added, it MUST filter by
  this contract (length / prefix) so it does not treat index pointers as
  task-entry records. The contract is documented at the top of the entry/index
  accessors in `DefaultCacheStore`.
- `CloudCacheStore` (`nf-cloudcache`): implement the same seam as a small object
  whose payload is the finalHash. Here separation is by **path** (a per-session
  `index/` subdir, see below), the natural namespacing for an object store —
  i.e. each backend uses the idiom native to its medium. This keeps the
  interface fully implemented across stores; no cloud coordination logic is
  added here.

**Session scoping (default / non-global cache).** This index is intentionally
**per-session**, preserving today's resume semantics — it is *not* a
cross-session/global index (that is the separate `nf-cloudcache-global`
feature, out of scope here). Scoping is enforced by *storage location*, not by
the key:

- `DefaultCacheStore` already isolates per session: its LevelDB lives at
  `<baseDir>/cache/<uniqueId>/db`, and `-resume` reuses the same `uniqueId`, so
  a prefixed key written there is inherently session-scoped.
- `CloudCacheStore` must place the index object under its per-session
  `dataPath` (`<basePath>/<uniqueId>`), i.e. at
  `<basePath>/<uniqueId>/index/<contentHash>` — **not** under a shared
  `<basePath>/index/...` prefix, which would leak pointers across sessions.

(The content hash already incorporates the session id via `TaskHasher`, so the
key is session-specific too; storage-location scoping makes the guarantee
explicit and matches how task entries are already partitioned.)

**2. `CacheDB`** (`nextflow.cache.CacheDB`) — thin delegators mirroring the
existing `getTaskEntry` / `putTaskAsync` pattern.

- `HashCode getSuccessfulHash(HashCode contentHash)` → `store.getSuccessfulHash(...)`.
- `void putSuccessfulHashAsync(HashCode contentHash, HashCode finalHash)` — enqueued
  on the **same `writer` Agent** used by `putTaskAsync`, so the index write is
  serialized *after* the corresponding entry write. This enforces the
  commit-marker discipline: a reader following the index never lands on an entry
  that is not yet durable.
- `writeTaskEntry0` writes a **fourth record slot** holding
  `task.contentHash?.asBytes()`, alongside the existing `[trace, ctx,
  refCount]`. This is a backward-compatible record extension: code that reads
  the record ignores extra slots, and entries written before this change simply
  lack slot 3.
- `removeTaskEntry(finalHash)`, when the ref-count reaches zero and the entry is
  deleted, recovers the content hash from slot 3 and removes the index pointer
  **iff** `getSuccessfulHash(contentHash) == finalHash` (the pointer targets the
  deleted entry). Old 3-slot entries → no content hash → the pointer is left to
  self-heal on the next resume.

**3. `TaskRun`** (`nextflow.processor.TaskRun`).

- New field `HashCode contentHash`, set from the `hash` argument inside
  `checkCachedOrLaunchTask`. Required because `task.hash` only ever holds the
  finalHash, and the index write at task-completion time needs the key.
- This field **replaces** PR #6903's `previousTryCount`; that field and its
  associated increment logic are not introduced.

### Read path (resume)

A minimal fast-path is added at the top of `checkCachedOrLaunchTask`, before the
existing loop — the content-hash capture plus a single delegation to a
`resumeFromSuccessfulHashIndex` helper. The loop body is left **byte-for-byte unchanged**;
all new lookup logic lives in the helper, keeping the method's complexity flat:

```
// Capture the content hash ONLY on the first attempt. Retries re-enter this
// method via resumeOrDie -> checkCachedOrLaunchTask(taskCopy, taskCopy.hash,
// false), where the `hash` argument is the previous attempt's (chained) final
// hash, not the content hash. makeCopy()/clone() carry contentHash forward, so
// the null guard preserves the original key across the whole retry chain.
if( task.contentHash == null )
    task.contentHash = hash

// fast resume path; otherwise fall back to the scan below
if( shouldTryCache && resumeFromSuccessfulHashIndex(task) )
    return

// existing iterate-and-bump loop, unchanged
```

with the helper:

```
private boolean resumeFromSuccessfulHashIndex( TaskRun task ) {
    try {
        final indexed = session.cache.getSuccessfulHash(task.contentHash)
        if( !indexed )
            return false
        final entry = session.cache.getTaskEntry(indexed, this)
        final resumeDir = entry ? FileHelper.asPath(entry.trace.getWorkDir()) : null
        return entry && entry.trace.isCompleted() && resumeDir?.exists()
                && checkCachedOutput(task.clone(), resumeDir, indexed, entry)
    }
    catch( Throwable t ) {
        log.trace("Hash-index lookup failed -- falling back to scan", t)
        return false   // stale/invalid pointer or error -> scan; self-heals on success
    }
}
```

**Why the hash sequence chains.** On the first call, `hash` is the content
hash and the loop computes `H(hash + tries)`, reassigning `hash` in place each
iteration — so attempt 2 is `H(H(contentHash + 1) + 2)`, etc. The resume walk
reproduces the same chain because it reassigns `hash` identically. The index
therefore stores `contentHash → <final chained hash of the successful
attempt>`, and the fast-path jumps straight to it without replaying the chain.

### Write path

The index write is centralized in `Session`, alongside the existing entry
writes, so it covers every success path without touching the
`checkCachedOrLaunchTask` loop body:

- **`notifyTaskComplete`** (a task finished executing): after
  `cache.putTaskAsync(handler, trace)`, if `trace?.isCompleted()` and
  `handler.task.contentHash` is set, enqueue
  `cache.putSuccessfulHashAsync(handler.task.contentHash, handler.task.hash)`.
- **`notifyTaskCached`** (a task resumed from cache — fired by
  `checkCachedOutput` for both the fast-path and the scan-path): after the
  existing `cache.cacheTaskAsync(handler)`, if `trace` is present and
  `handler.task.contentHash` is set, enqueue the same
  `putSuccessfulHashAsync(...)`. This rewrites the pointer on every resume, which
  both self-heals a stale pointer and upgrades a scan-path resume to a
  fast-path resume next time.

Because both writes go through the `writer` Agent *after* the corresponding
entry write, the commit-marker discipline holds. Writes are last-writer-wins:
the process is single-writer, and all successful executions for a given
contentHash are content-equivalent, so the write is idempotent.

### Edge cases

| Case | Behaviour |
|------|-----------|
| Stale pointer (entry GC'd, workdir gone, output check fails) | Fast-path validation fails → fall through to loop → success rewrites the pointer. Self-heals in one lookup. |
| Pointer targets an incomplete entry | Cannot occur under commit-marker discipline; defensive `isCompleted()` check treats it as stale regardless. |
| Two tasks, same contentHash, within one run | `tries++` + `lockManager` allocate distinct workdirs as today; index write is idempotent. No new race. |
| Retry within the same run | Index is written only on success; a failed attempt writes nothing. The relaunch succeeds at a fresh `tries` and writes the pointer. No `previousTryCount` needed. |
| `getSuccessfulHash` throws / corrupt pointer | Caught; logged at trace; fall through to loop. The index is strictly an optimization, never fatal. |
| `removeTaskEntry` deletes the **successful** entry (ref-count → 0) | The content hash is recovered from record slot 3; since `getSuccessfulHash(contentHash) == finalHash`, the pointer is removed too. |
| `removeTaskEntry` deletes a **failed-attempt** entry sharing the same content hash | `getSuccessfulHash(contentHash)` points at the *successful* finalHash (≠ the deleted hash), so the pointer is **preserved**. |
| `removeTaskEntry` deletes an old 3-slot entry (no content hash) | No content hash to recover → pointer untouched; self-heals on next resume. |
| `nextflow clean <runName>` of a single-run session | Every entry is removed and then `drop()` wipes the whole per-session DB, index included. |

### Backward compatibility

- Old caches have no index keys → `getSuccessfulHash` returns `null` → the existing
  loop runs exactly as today. Such caches retain the original #6884 limitation
  until a fresh run repopulates the index, after which resume is correct.
- No on-disk format version bump: new namespaced index keys are simply absent in
  old LevelDB databases, and the new fourth record slot is a backward-compatible
  extension (readers ignore extra slots; old 3-slot entries lack it).
- Nothing to migrate from #6903, since it is not merged.
- An entry written before this change carries no content hash, so deleting it
  cannot clean its index pointer; the orphaned pointer self-heals on the next
  resume. A full cache `drop()` removes any remaining pointers.

## Testing

- **`DefaultCacheStoreTest`** — round-trip `putSuccessfulHash` / `getSuccessfulHash` /
  `deleteSuccessfulHash`; absent key returns `null`; prefixed index keys never
  collide with entry keys; `drop()` removes index keys.
- **`CacheDBTest`** — delegators; async ordering (entry committed before its
  index pointer); `removeTaskEntry` removes the pointer for the successful
  entry it targets, **preserves** the pointer when a failed-attempt entry
  sharing the same content hash is removed.
- **`TaskProcessorTest`** (primary behavioural coverage):
  - success after N failed attempts resolves on resume in a single index
    lookup, with no scan;
  - the exact #6884 reproduction (fail → retry → succeed, kill, resume, kill,
    resume) resumes the completed task instead of re-executing it;
  - stale pointer (target workdir removed) falls back, re-resumes or
    re-executes, and rewrites the pointer;
  - index miss on a legacy entry falls through to the loop unchanged.
- **`CloudCacheStoreTest`** — `index/<hash>` object round-trip against a mocked
  store.
- **Integration** — a `tests/` case reproducing #6884 end-to-end.

## Files affected

- `modules/nextflow/src/main/groovy/nextflow/cache/CacheStore.groovy`
- `modules/nextflow/src/main/groovy/nextflow/cache/DefaultCacheStore.groovy`
- `modules/nextflow/src/main/groovy/nextflow/cache/CacheDB.groovy`
- `modules/nextflow/src/main/groovy/nextflow/processor/TaskProcessor.groovy`
- `modules/nextflow/src/main/groovy/nextflow/processor/TaskRun.groovy`
- `plugins/nf-cloudcache/src/main/nextflow/cache/CloudCacheStore.groovy`
- Tests as listed above.