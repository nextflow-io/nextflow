# Spec-driven versioned task hasher

- Author: Jorge Ejarque
- Status: draft — design agreed, prototype not yet built
- Date: 2026-09-16
- Related: `adr/20260901-task-hash-key-fields.md`, nextflow-io/nextflow#6927, seqeralabs/nextflow#47

## Summary

Replace the family of `TaskHasher` subclasses with **one interpreter plus data**. A `TaskHashSpec`
declares an ordered list of named key bindings, a set of encoding rules and a hash function;
`BaseTaskHasher` executes any spec. Hash versions, the global-cache hasher and any future variant
become spec values rather than classes.

Two capabilities fall out, and they are the reason for doing it:

1. **Per-key digests** — `explain()` returns `key name → digest` for a task, so two runs can be
   diffed by key name even across hash versions. This is what the Seqera Platform cache-miss feature
   in `adr/20260901-task-hash-key-fields.md` needs.
2. **A governed key vocabulary** — one enum of key names, shared by the hasher, the lineage record
   mapping and `-dump-hashes`, with an exhaustiveness test so a new key cannot be added silently.

## Problem

Three independent dimensions of variation exist in code today, and subclassing multiplies them:

| | nextflow-io#6927 | seqeralabs#47 |
| --- | --- | --- |
| key **set** | V3 adds project `bin/` entries | drops `sessionId` and `processName` |
| key **encoding** | `withOrderIndependentMaps`, `withCacheFunnelFirst` | — |
| key **extraction** | — | file inputs by content identity |

`GlobalTaskHasher extends TaskHasher` and re-implements the whole of `compute()` — roughly 60 lines
duplicated to express two differences. A key added to core is silently absent from the global hasher.

Nothing identifies which key set produced a given hash, so a consumer cannot tell whether two hashes
are even comparable.

## Decisions

| # | Decision | Rationale |
| --- | --- | --- |
| D1 | **Byte-exact reproduction is the acceptance test.** Any spec representing an existing hasher must produce an identical hash for the same task. | A task executed under an older release, or under the global cache, must still hit the cache when re-run under the new mechanism. Rules out otherwise-desirable cleanups inside existing specs. |
| D2 | Reproduce **all four historical eras** of master's hash (H1–H4) plus #47's global hasher. | Exercises all three variation dimensions against real code. Note this supersedes the original framing of "reproduce #6927's V1/V2/V3": see [Finding](#finding-6927s-versions-do-not-map-onto-history) — two of those three correspond to no release that ever shipped. |
| D3 | The hasher must emit **per-key digests plus the spec** (id + declared key names), not just the final hash. | Diffing by key name works even where the lineage record stores a lossy form of the value (`codeChecksum` vs `source`). |
| D4 | Version identity = **declared name + derived fingerprint**, with a test asserting they move together. | Readable in a UI; impossible to forget; stable under refactoring because the canonical form excludes implementation detail. |
| D5 | **Central closed `HashKey` enum**; specs select, order and bind. | Diff alignment is by key name, so names are load-bearing data. Free strings would let two specs drift in naming and would move a historical spec's fingerprint on an innocuous rename. |
| D6 | Validation is **end-to-end via `-resume`**, using a dedicated pipeline. | A cache hit on resume *is* byte-exact hash equality, proven through `CacheDB` and `TaskProcessor` rather than a unit assertion. |
| D7 | Pipeline covers container and conda; **spack and the `module` directive** are unit-oracle only. | Those need HPC/solver infrastructure; they are simple pass-through values. |
| D8 | **Opt-in only**, through master's merged `TaskHasherFactory`. No version requested → the factory abstains → the stock `TaskHasher` runs, untouched. | Removes the design's largest risk outright: if the default path cannot change, it cannot invalidate anyone's cache. Also retires the `TaskHashSpecFactory` SPI and every `TaskProcessor` change this spec previously required. Added 2026-09-23 after master landed the hook. |
| D9 | **Specs are JSON resources**, resolved against a `ContributorRegistry`; the loader validates strictly and never falls back. | Composition is data, extraction is code. A version that only reorders or re-includes keys ships as a file, and the file doubles as the ADR's hash-relevant field list. Strict validation because a mis-resolved spec produces cache keys nobody can account for. Added 2026-09-23. |

## Architecture

### Types

```groovy
// nf-commons — the encoding half of the hash version
class EncodingRules {
    boolean orderIndependentMaps      // Map → values() (false) | hashUnorderedCollection(entrySet()) (true)
    boolean cacheFunnelFirst          // CacheFunnel dispatched before (true) | after (false) Map/Set
}
```

```groovy
// nextflow.processor.hash — the key half
enum HashKey {
    SESSION_ID, PROCESS_NAME, TASK_SOURCE, CONTAINER, INPUTS, EVAL_OUTPUTS,
    SCRIPT_VARS, BIN_ENTRIES, MODULE_BUNDLE, ENV_MODULES, CONDA, SPACK, STUB_MARKER
}

interface Contributor {
    List<Object> emit(HashContext ctx)     // 0..N values in order; [] means "contributes nothing"
}

class KeyBinding { HashKey key; Contributor contributor }

class TaskHashSpec {
    String id                              // 'std/v4', 'global/v1'
    List<KeyBinding> bindings              // ordered
    EncodingRules encoding
    HashFunction function
    String fingerprint()                   // over canonical form: id + ordered key names + encoding + function
}

class HashContext { TaskRun task; TaskProcessor processor; Session session; Map services }

class BaseTaskHasher implements TaskHasher {
    HashCode compute()                     // fold every emitted value, in order
    Map<HashKey,HashCode> explain()        // per-binding digest — additive, never feeds compute()
}
```

### Contributor versus EncodingRules

A `Contributor` chooses **which object** enters the hash. `EncodingRules` decide **how any object
becomes bytes**. Test: feed in the identical object — if the bytes differ, it is encoding.

The split is forced, not stylistic. A contributor can only choose the *top-level* object it emits;
`HashBuilder` traverses recursively, so its rules also govern objects nested inside — a `Map` buried
in an input value, a collection inside a `CacheFunnel`. That is exactly what #6679 changed, and no
contributor can express it.

The borderline case is #47's `contributeInput`, which recurses into nested collections replacing file
values with identity tokens. It *substitutes different objects*, so it is a contributor despite
traversing.

Rule of thumb when writing specs: if today's code makes the choice at the call site, it is a
contributor; if `HashBuilder` makes it during traversal, it is encoding.

Scope: `EncodingRules` is per-spec, not per-key. Per-key encoding reproduces nothing that exists.

### Byte-exactness hazards

All of them live in the contributors:

- **Emission arity.** `INPUTS` emits two values per input (name, then value); `BIN_ENTRIES` and
  `ENV_MODULES` use `addAll` (N values); `EVAL_OUTPUTS` emits the literal `"eval_outputs"` marker
  *plus* the payload.
- **Conditional emission must emit nothing, not null.** Every optional key appends zero elements when
  absent. Returning `[null]` changes the list length and therefore the hash.
- **`SPACK` emits `[spack]` or `[spack, arch]`.** `arch` sits inside the `spack` conditional in the
  current code, so it is not an independent key. (Note for the ADR: the drift table lists
  `architecture` as a record field, but it only contributes when spack is set.)
- **`MODULE_BUNDLE` is doubly conditional** — `session.enableModuleBinaries()` *and*
  `bundle.hasEntries()`.
- **`HashMode` is resolved at compute time, not held on the spec.** `compute()` takes no arguments
  (matching the existing `TaskHasher` interface) and reads
  `task.processor.getConfig().getHashMode()` exactly as the current code does. The mode is a
  per-process runtime choice (`cache 'deep'`), orthogonal to the spec — two runs can share a spec and
  hash differently under `DEEP`.

## The spec matrix

### Verified history of master's key list

Traced from `TaskHasher.groovy` and `HashBuilder.java`:

| Era | Window | eval form | `bin/` entries | module bundle | encoding |
| --- | --- | --- | --- | --- | --- |
| H1 | … → 2026-03-09 | derived string | yes | no | pre-record-types |
| H2 | 2026-03-09 → 2026-07-17 | derived string | yes | no | record-types |
| H3 | 2026-07-17 → 2026-09-03 | derived string | yes | yes | record-types |
| H4 | 2026-09-03 → today | **raw map** | yes | yes | record-types |

Boundary commits:

- `d54ff29af` #6679 "Record types" (2026-03-09) — changed `HashBuilder` in exactly two ways: `Map`
  from `values()` to `hashUnorderedCollection(entrySet())`, and `CacheFunnel` dispatched *before*
  `Map`/`Set`. These are precisely #6927's two flags.
- `029e52eef` #6914 "Fix cache invalidation for module binaries" (2026-07-17) — added the
  module-bundle key.
- `9e7a492a3` #7575 "Make eval output names stable across runs" (2026-09-03) — replaced
  `computeEvalOutputCommands(outEvals)` with `outEvals` and deleted the method.

`bin/` entries have been hashed since `516294e7e` (2018). **There has never been an era without
them.**

### Finding: #6927's versions do not map onto history

| #6927 | `bin/` | bundle | encoding | eval | matches |
| --- | --- | --- | --- | --- | --- |
| V1 | no | yes | pre-record-types | derived | **nothing** |
| V2 | no | yes | record-types | derived | **nothing** |
| V3 | yes | yes | record-types | derived | H3 |
| — | | | | | H4 (master today) has no version |

Only V3 reproduces a real historical hash. V1 and V2 are points in a flag space: they omit `bin/`
entries, which no release ever omitted, while including the module bundle, which no pre-July release
ever had. V1's stated purpose — "reproduce pre-record-types hashes for cache backward compatibility"
— is not achieved.

Separately, #6927 predates #7575, so **its default V3 reverts that change**: merging as-is silently
invalidates the cache for every pipeline using eval outputs. Both points should be raised on the PR.

### Specs to implement

Anchored to history rather than to flags.

| `HashKey` | `std/v1`=H1 | `std/v2`=H2 | `std/v3`=H3 | `std/v4`=H4 (master) | `global/v1` (#47) |
| --- | --- | --- | --- | --- | --- |
| `SESSION_ID` | yes | yes | yes | yes | **no** |
| `PROCESS_NAME` | yes | yes | yes | yes | **no** |
| `TASK_SOURCE` | yes | yes | yes | yes | yes |
| `CONTAINER` | yes | yes | yes | yes | yes |
| `INPUTS` | raw value | raw | raw | raw | **file identity** |
| `EVAL_OUTPUTS` | derived string | derived string | derived string | **raw map** | raw map |
| `SCRIPT_VARS` | yes | yes | yes | yes | yes |
| `BIN_ENTRIES` | yes | yes | yes | yes | yes |
| `MODULE_BUNDLE` | **no** | **no** | yes | yes | yes |
| `ENV_MODULES` | yes | yes | yes | yes | yes |
| `CONDA` | yes | yes | yes | yes | yes |
| `SPACK` | yes | yes | yes | yes | yes |
| `STUB_MARKER` | yes | yes | yes | yes | yes |
| `orderIndependentMaps` | false | **true** | true | true | true |
| `cacheFunnelFirst` | false | **true** | true | true | true |

Only two cells differ between adjacent `std` versions: the specs are diffs, not copies.

## Integration and selection

### Call site — now supplied by master

`TaskHasherFactory` and `TaskProcessor.createTaskHasher()` **landed on master** while this prototype
was being built (the branch was rebased onto it on 2026-09-23). Master asks plugin factories in
priority order, takes the first non-null hasher, and otherwise uses the stock `TaskHasher` — and it
caches the factory list in a `@PackageScope volatile List<TaskHasherFactory>` field on
`TaskProcessor`, arriving independently at the same publication fix this prototype had made.

Consequence: the `TaskHashSpecFactory` SPI and the `TaskProcessor` wiring this spec previously
described are **superseded and were dropped**. Nothing in this design touches `TaskProcessor` now.

### Seam: opt in through master's hook

`SpecTaskHasher extends TaskHasher` and overrides `compute()`. `SpecTaskHasherFactory implements
TaskHasherFactory` returns one **only when a version is requested**, and `null` otherwise:

```groovy
TaskHasher create(TaskRun task) {
    final spec = TaskHashSpecResolver.requestedSpec()   // null unless NXF_TASK_HASH_VER is set
    return spec == null ? null : new SpecTaskHasher(task, spec)
}
```

It is registered in core's `META-INF/extensions.idx` beside `DefaultCacheFactory` and
`DefaultTaskTipProvider`, so core needs no plumbing of its own.

This replaces the earlier "plugins contribute specs, not hashers" seam, and is better on the point
that mattered most:

- **The default path is untouched.** With no version requested the factory abstains and the stock
  `TaskHasher` runs. Opting in is the only way to change a cache key, so the design carries no risk
  of invalidating anyone's cache — which was the single largest objection to it.
- **`std/v4` becomes a permanent equivalence oracle** rather than the production path. The test
  asserting `SpecTaskHasher(STD_V4) == TaskHasher` on fixtures is now the drift guard (A2) the ADR
  asks for, with a standing reason to exist.
- What the old seam bought — that a plugin cannot bypass `explain()` — is given up. A plugin can
  still subclass `TaskHasher` directly, as `GlobalTaskHasher` does. That is master's call, not this
  design's, and re-litigating it is not worth the coupling.

### Specs are JSON, not code

A spec is composition — which keys, in what order, under which encoding and hash function — and only
the *extraction* is code. So a spec is a JSON resource:

```json
{ "id": "std/v3",
  "function": "murmur3_128",
  "encoding": { "orderIndependentMaps": true, "cacheFunnelFirst": true },
  "keys": [ { "key": "SESSION_ID", "contributor": "sessionId" },
            { "key": "EVAL_OUTPUTS", "contributor": "evalOutputs.derivedString" } ] }
```

`ContributorRegistry` resolves `contributor` names; plugins register their own so a plugin spec can
name contributors core does not know about. `TaskHashSpecLoader` validates strictly — unknown key,
unknown contributor, duplicate binding, unsupported function, malformed document all fail at load,
because a spec silently resolving to something else produces cache keys nobody can account for.

Two consequences worth stating:

- **A version that only adds, removes or reorders keys, or flips an encoding flag, needs no Nextflow
  release.** #6914 and #6679 would each have been one new file. A version needing a *new contributor*
  still needs code, as #7575 did — JSON covers composition, not derivation.
- **Diffing two spec files answers "what changed between these hash versions"**, which is the same
  artifact `adr/20260901-task-hash-key-fields.md` calls the hash-relevant field list. One artifact,
  two uses — and it makes the ADR's axis 2 concrete: shipped in the jar is B2-shaped, a separate
  "hash version repository" is B4-shaped with B4's unresolvable-for-unknown-builds failure. Shipping
  them in the jar *and* mirroring the repo gets offline resolution plus a correction layer.

### Spec selection

`NXF_TASK_HASH_VER` names a version; absent, no spec is used at all. An unknown id is an error, never
a silent fallback.

### Placement

| Module | Contents |
| --- | --- |
| `nf-commons` | `EncodingRules`; the two `HashBuilder` flags, defaults preserving master |
| `nextflow.processor.hash` | `HashKey`, `Contributor`, `KeyBinding`, `HashContext`, `TaskHashSpec`, `TaskHashSpecLoader`, `ContributorRegistry`, `Contributors`, `SpecTaskHasher`, `SpecTaskHasherFactory`, `TaskHashSpecResolver`, `StdSpecs` |
| `nextflow` resources | `nextflow/processor/hash/std-v{1..4}.json`, plus the factory line in `META-INF/extensions.idx` |
| `nf-cloudcache-global` | the `global/v1` spec |

Reproducing `std/v1` and `std/v2` requires the `HashBuilder` encoding flags, which exist only in
#6927. The prototype ports that change itself, defaulting to current behaviour so master is
unaffected — deliberate overlap that demonstrates how the two reconcile.

### Immediate payoff, no record change

Wire `explain()` into `-dump-hashes json`, which today emits `hash`/`type`/`value` with no names.
Named entries make the primary cache-miss debugging tool readable, and exercise the diff path on
every run without touching the lineage model — which stays the ADR's undecided question.

### Default-unchanged guarantee

Since D8 this is structural rather than demonstrated: with no version requested, no spec is
constructed and `TaskProcessor` runs the stock `TaskHasher`. The byte-exactness evidence below is
still what makes `std/v4` trustworthy *as an oracle*, but it is no longer what protects the default
path.

#### Historical note

No env var and no plugin → `std/v4` → byte-identical to master. This is test 1 of the runbook.

## Validation

### Baselines

The era map lands on released edge builds, so three of five baselines need no compilation:

| Spec | Baseline | Why |
| --- | --- | --- |
| `std/v1` | `NXF_VER=26.02.0-edge` (2026-02-28) | before #6679 |
| `std/v2` | `NXF_VER=26.07.0-edge` (2026-07-15) | after #6679, before #6914 |
| `std/v3` | `NXF_VER=26.08.0-edge` (2026-08-20) | after #6914, before #7575 |
| `std/v4` | master build | current |
| `global/v1` | #47 build | — |

### Pipeline

`specs/260916-task-hash-spec/pipeline/` — one process per key, so a miss localises immediately.

| Process | Exercises | Note |
| --- | --- | --- |
| `P_BASIC` | `TASK_SOURCE`, `INPUTS` (value + file), `SCRIPT_VARS` (global var + `task.ext`) | |
| `P_MAP_INPUT` | encoding | **H1↔H2 discriminator** — needs a genuine `Map` input |
| `P_EVAL` | `EVAL_OUTPUTS` | **H3↔H4 discriminator** |
| `P_MODULE_BUNDLE` | `MODULE_BUNDLE` | **H2↔H3 discriminator** — needs module binaries enabled |
| `P_BIN` | `BIN_ENTRIES` | calls a `bin/` script |
| `P_CONTAINER` | `CONTAINER` | image pinned by digest |
| `P_CONDA` | `CONDA` | env pinned |
| `P_STUB` | `STUB_MARKER` | run with `-stub-run` |
| `P_PRODUCER` → `P_CONSUMER`, plus an external file | file identity | `global/v1` only; includes a nested list input for `contributeInput` recursion |

`SESSION_ID` and `PROCESS_NAME` are exercised implicitly — resume restores the session id, and the
negative controls cover the rest. `ENV_MODULES` and `SPACK` are unit-oracle only.

### Three assertions, not one

1. **Positive** — baseline run, then prototype with `NXF_TASK_HASH_VER=<spec> -resume`: every task
   reports `CACHED`.
2. **Negative controls** — change a script body, rename a process, change a value input: those tasks
   must re-execute. Without these, a run where everything caches for unrelated reasons reads as a
   pass.
3. **Cross-spec discrimination** — resume an H3 baseline with `std/v4` and everything caches
   **except `P_EVAL`**; resume an H2 baseline with `std/v3` and only `P_MODULE_BUNDLE` misses. This
   proves the specs are genuinely distinct in exactly the predicted place, not merely
   self-consistent.

### Confounds

- **Cache-DB compatibility.** Resuming a February build's LevelDB/Kryo cache with a master build may
  fail for reasons unrelated to hashing. Mitigation: an independent second track — run both builds
  with `-dump-hashes json` and compare the reported task hashes directly. It isolates hashing from
  cache plumbing, and per-key entries localise any mismatch. Use it whenever resume fails.
- **Environment drift.** Container digests and resolved conda environments must be pinned, or
  `CONTAINER` and `CONDA` differ for environmental reasons.

### Unit oracles

For `ENV_MODULES` and `SPACK`, and for any key the pipeline cannot reach: assert
`BaseTaskHasher(spec).compute(task) == referenceImplementation(task)` over a fixture matrix, with the
reference copied read-only into test sources.

## Risks and open items

- **The era map derived from `TaskHasher` and `HashBuilder` alone is incomplete — now DEMONSTRATED,
  not hypothetical.** Validation against genuine releases (2026-09-22) found a fourth boundary:
  `785e801ad` #7165 "Fix strict parser to include params refs in task hash" (2026-05-21), which lives
  in the parser/nf-lang layer. It makes `getTaskGlobalVars()` fold referenced `params.*` into
  `SCRIPT_VARS`. `std/v2` and `std/v3` reproduce their releases 8/8 because both post-date it;
  `std/v1` misses exactly the one process that reads `params`, because `26.02.0-edge` predates it.

- **A spec cannot express how an era DERIVED a hashed value — only which keys are hashed and how they
  are encoded.** This is the sharpest limit the validation found, and it was invisible until the
  prototype was run against a real old release. Contributors call *today's* helpers: `SCRIPT_VARS`
  delegates to `getTaskGlobalVars()`, `CONTAINER` to `getContainerFingerprint()`, `CONDA` to
  `getCondaEnv()`. If any of those changes what it returns, every spec's output moves with it and no
  `TaskHashSpec` can hold the old behaviour. Reproducing an era across such a change would require
  era-appropriate *helpers*, not just an era-appropriate key list — a strictly larger design.
  Practical consequence: specs reproduce eras reliably only back to the most recent input-derivation
  change, which bounds how far back the cache-miss explanation can honestly reach.
- **`std/v1` and `std/v2` composition is a hypothesis** until the resume test confirms it.
- **Fingerprint canonical form must be frozen carefully.** It has to include extraction identity —
  `global/v1` with `SampleFileIdentity` and with `ProvenanceFileIdentity` must not share a
  fingerprint — while excluding implementation detail, so historical ids stay stable.
- **Relationship to #6927 is unresolved.** If it lands, the two must be reconciled rather than
  coexisting as separate identifiers. The findings above should be raised on that PR first.
- **Lineage storage is out of scope.** What the record carries remains
  `adr/20260901-task-hash-key-fields.md`'s decision; this spec only guarantees the hasher can emit
  it.
- **A document-based hash was considered and rejected.** Building a canonical serialised document and
  hashing it (Nix-derivation style) has the most explanatory power, but today's hash is a streaming
  fold over heterogeneous objects, so a document hash cannot reproduce `std/*` or `global/*`. It
  remains available as a *future profile*, never as a reproduction.

## References

- `adr/20260901-task-hash-key-fields.md` — the lineage/Platform cache-miss ADR this enables
- nextflow-io/nextflow#6927 — versioned task hasher strategy (open since 2026-03)
- seqeralabs/nextflow#47 — global cache plugin, `TaskHasherFactory` SPI, `GlobalTaskHasher`
- nextflow-io/nextflow#6679 `d54ff29af`, #6914 `029e52eef`, #7575 `9e7a492a3` — the three boundaries
- `modules/nextflow/src/main/groovy/nextflow/processor/TaskHasher.groovy`
- `modules/nf-commons/src/main/nextflow/util/HashBuilder.java`
