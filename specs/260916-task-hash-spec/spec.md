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
| D2 | Reproduce **every historical era** of master's hash back to the supported floor, plus #47's global hasher. Superseded by D10: the sweep to 25.10.0 found seven eras, not four. | Exercises all three variation dimensions against real code. Note this supersedes the original framing of "reproduce #6927's V1/V2/V3": see [Finding](#finding-6927s-versions-do-not-map-onto-history) — two of those three correspond to no release that ever shipped. |
| D3 | The hasher must emit **per-key digests plus the spec** (id + declared key names), not just the final hash. | Diffing by key name works even where the lineage record stores a lossy form of the value (`codeChecksum` vs `source`). |
| D4 | Version identity = **declared name + derived fingerprint**, with a test asserting they move together. | Readable in a UI; impossible to forget; stable under refactoring because the canonical form excludes implementation detail. |
| D5 | **Central closed `HashKey` enum**; specs select, order and bind. | Diff alignment is by key name, so names are load-bearing data. Free strings would let two specs drift in naming and would move a historical spec's fingerprint on an innocuous rename. |
| D6 | Validation is **end-to-end via `-resume`**, using a dedicated pipeline. | A cache hit on resume *is* byte-exact hash equality, proven through `CacheDB` and `TaskProcessor` rather than a unit assertion. |
| D7 | Pipeline covers container and conda; **spack and the `module` directive** are unit-oracle only. | Those need HPC/solver infrastructure; they are simple pass-through values. |
| D8 | **Opt-in only**, through master's merged `TaskHasherFactory`. No version requested → the factory abstains → the stock `TaskHasher` runs, untouched. | Removes the design's largest risk outright: if the default path cannot change, it cannot invalidate anyone's cache. Also retires the `TaskHashSpecFactory` SPI and every `TaskProcessor` change this spec previously required. Added 2026-09-23 after master landed the hook. |
| D9 | **Specs are JSON resources**, resolved against a `ContributorRegistry`; the loader validates strictly and never falls back. | Composition is data, extraction is code. A version that only reorders or re-includes keys ships as a file, and the file doubles as the ADR's hash-relevant field list. Strict validation because a mis-resolved spec produces cache keys nobody can account for. Added 2026-09-23. |
| D10 | **25.10.0 is the oldest supported version**, giving seven THVs. | Lineage shipped in 25.04 and no record predates it, so there is nothing older to explain; 25.10.0 is a stable release inside that window. Added 2026-10-02 after sweeping 25.10.0→master, which found three boundaries the earlier four-version map had missed. |

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

### Verified history, with 25.10.0 as the floor

The floor is **25.10.0** (2025-10-22). Two reasons: lineage shipped in 25.04 and is still marked
experimental, so no lineage record predates April 2025; and 25.10.0 is a stable release inside that
window, which makes it a defensible oldest-supported point.

Swept from 25.10.0 to master across `TaskHasher`, `HashBuilder`, `CacheHelper`, the value-derivation
helpers, and the parser. Six boundaries, so **seven versions**:

| Ver | First release | Boundary that ends it | Hash differs only when |
| --- | --- | --- | --- |
| `std/v1` | 25.10.0 | `1ca327c80` #6605 | a hashed file sits under the asset root but not under `baseDir` — running from a repo subdirectory or with `-main-script` |
| `std/v2` | 25.11.0-edge | `b5278c75a` #6696 | the script references `task.ext.*` |
| `std/v3` | 26.01.1-edge | `d54ff29af` #6679 | a hashed value is a `Map` |
| `std/v4` | 26.03.0-edge | `785e801ad` #7165 | the script references `params.*` in a process body |
| `std/v5` | 26.04.2 | `029e52eef` #6914 | module binaries enabled **and** the bundle has entries |
| `std/v6` | 26.08.0-edge | `9e7a492a3` #7575 | the process declares `eval` outputs |
| `std/v7` | 26.09.0-edge | — (current) | — |

Most pipelines cross most boundaries with no hash change. A pipeline with no evals, no map inputs,
no `params.*` or `task.ext.*` in a process body, no module bundle, run from the repo root, hashes
identically from 25.10.0 to today.

Two boundaries are **not** key-set or encoding changes. #6696 and #7165 change what a helper
*returns* — the value, not the key list. They are still cache boundaries, and both are reproducible
because each one *adds entries carrying an identifying prefix*, so an older spec recovers the old
value by filtering: drop `task.ext.*` for `std/v1`, drop `params.*` for `std/v1`–`std/v4`.

### Measured and rejected as a boundary: #6789

`66b743836` #6789 "Fix different task hash with v2 parser" (2026-02-03) looked like an
unreproducible boundary, because it changes extracted source text rather than adding entries. It is
not a boundary at all, and the measurement is worth recording.

Comparing the **task-source entry** across runs (not the whole task hash — that embeds
`session.uniqueId`, which is new per run, so whole-hash comparison across runs always differs and
proves nothing):

| Comparison | Result |
| --- | --- |
| 25.10.4, parser v1 vs v2 | all 8 processes differ — the bug #6789 fixed |
| current master, parser v1 vs v2 | all 8 identical — the fix holds |
| 25.10.4 v1 vs current master v1 | all 8 identical — v1 source text never moved |

The v1 parser was the default until `295f17307` #6643 (2026-01-16), and #6789 landed 2026-02-03.
Both ship first in the same release, 26.01.1-edge, so **no released default-parser build ever
produced the broken v2 source text**. The boundary only existed for someone who opted into v2
in that window.

### What is implemented today

Four of the seven exist. The current ids are off by the renumbering above:

| New id | Current id | State |
| --- | --- | --- |
| `std/v1` | — | not implemented; needs an asset-detection encoding flag and a `task.ext.*` filter |
| `std/v2` | — | not implemented; needs a `task.ext.*` filter |
| `std/v3` | `std/v1` | implemented, validated 8/8 against 26.02.0-edge |
| `std/v4` | — | not implemented; needs a `params.*` filter |
| `std/v5` | `std/v2` | implemented, validated 8/8 against 26.04.6 |
| `std/v6` | `std/v3` | implemented, validated 8/8 against 26.08.0-edge |
| `std/v7` | `std/v4` | implemented, equals current master |

Renumbering changes every spec id, and therefore every fingerprint. Nothing persists a fingerprint
yet, so the cost is editing four JSON files and their tests.

The `7/8` result previously recorded for `std/v1` against 26.02.0-edge was this gap: that release
predates #7165, and the spec had no `params.*` filter. Under the new numbering that spec is
`std/v3`, and the missing filter belongs to `std/v4` and below.

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

Seven `std` specs plus the plugin's `global/v1`. `SESSION_ID`, `PROCESS_NAME`, `TASK_SOURCE`,
`CONTAINER`, `SCRIPT_VARS`, `BIN_ENTRIES`, `ENV_MODULES`, `CONDA`, `SPACK` and `STUB_MARKER` are
present in every `std` spec and are omitted from the table — only the cells that move are shown.

| | `v1` | `v2` | `v3` | `v4` | `v5` | `v6` | `v7` | `global/v1` |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `SESSION_ID` / `PROCESS_NAME` | yes | yes | yes | yes | yes | yes | yes | **no** |
| `INPUTS` | raw | raw | raw | raw | raw | raw | raw | **file identity** |
| `EVAL_OUTPUTS` | derived | derived | derived | derived | derived | derived | **raw map** | raw map |
| `MODULE_BUNDLE` | no | no | no | no | no | **yes** | yes | yes |
| `orderIndependentMaps` / `cacheFunnelFirst` | false | false | false | **true** | true | true | true | true |
| `SCRIPT_VARS` filter | `task.ext.*` + `params.*` | `task.ext.*` + `params.*` | `params.*` | `params.*` | none | none | none | none |
| asset detection | narrow | wide | wide | wide | wide | wide | wide | wide |

Adjacent specs differ in exactly one row. That is the property the cross-spec validation checks: a
run resumed one version back must miss only the process that exercises that row.

Three new mechanisms are needed beyond what is built:

- **`task.ext.*` filter** on `SCRIPT_VARS`, for `v1`–`v2` — #6696 ends `v2`, so neither collects
  `task.ext.*`. Same shape as the `params.*` filter.
- **`params.*` filter** on `SCRIPT_VARS`, for `v1`–`v4`.
- **Asset-detection flag** in `EncodingRules`, for `v1` only — the narrow `isAssetFile` that checked
  `baseDir` alone, before #6605 widened it to the asset root.

The first two are contributors. The third is a third `HashBuilder` flag, which also means
`EncodingRules.canonicalForm()` gains a term — safe to change now, since no fingerprint is persisted
anywhere yet.

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

`VersionedTaskHasher extends TaskHasher` and overrides `compute()`. `VersionedTaskHasherFactory implements
TaskHasherFactory` returns one **only when a version is requested**, and `null` otherwise:

```groovy
TaskHasher create(TaskRun task) {
    final spec = TaskHashSpecResolver.requestedSpec()   // null unless NXF_TASK_HASH_VER is set
    return spec == null ? null : new VersionedTaskHasher(task, spec)
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
  asserting `VersionedTaskHasher(STD_V4) == TaskHasher` on fixtures is now the drift guard (A2) the ADR
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
| `nextflow.processor.hash` | `HashKey`, `Contributor`, `KeyBinding`, `HashContext`, `TaskHashSpec`, `TaskHashSpecLoader`, `ContributorRegistry`, `Contributors`, `VersionedTaskHasher`, `VersionedTaskHasherFactory`, `TaskHashSpecResolver`, `StdSpecs` |
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

One genuine release per version. **Check the build number first** — `NXF_VER=<v> nextflow -version`
must show a real build; `build 0` means a locally installed snapshot, not the release, and such a
build will produce a wrong result (see Risks).

| Spec | Baseline | Build | Status |
| --- | --- | --- | --- |
| `std/v1` | 25.10.0 … 25.10.2 | — | not implemented, not run |
| `std/v2` | 25.12.0-edge | `build 0` locally — needs a clean download | not implemented, not run |
| `std/v3` | 26.02.0-edge | 11371 | **8/8 cached** (run as the old `std/v1`) |
| `std/v4` | 26.03.x … 26.04.1 | — | not implemented, not run |
| `std/v5` | 26.04.6 | 12646 | **8/8 cached** (run as the old `std/v2`) |
| `std/v6` | 26.08.0-edge | 13213 | **8/8 cached** (run as the old `std/v3`) |
| `std/v7` | master, or 26.09.1-edge | — | matches master; not yet run against a real release |
| `global/v1` | seqeralabs#47 build | — | not implemented |

Two gaps worth naming: `std/v7` has never been checked against a genuine release — 26.09.1-edge
would do it — and `std/v2`'s nearest baseline, 25.12.0-edge, is a `build 0` copy locally.

### Pipeline

`specs/260916-task-hash-spec/pipeline/` — one process per key, so a miss localises immediately.

| Process | Exercises | Note |
| --- | --- | --- |
| `P_BASIC` | `TASK_SOURCE`, `INPUTS` (value + file), `SCRIPT_VARS` (global var + `task.ext`) | |
| `P_MAP_INPUT` | encoding | **`v3`↔`v4` discriminator** — needs a genuine `Map` input |
| `P_EVAL` | `EVAL_OUTPUTS` | **`v6`↔`v7` discriminator** |
| `P_MODULE_BUNDLE` | `MODULE_BUNDLE` | **`v5`↔`v6` discriminator** — needs module binaries enabled |
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

- **A derivation change is reversible only when it ADDS identifiable entries.** #6696 and #7165 both
  add variable refs carrying a prefix (`task.ext.`, `params.`), so an older spec recovers the old
  value by filtering. #6789 changed the *content* of the extracted source text instead, which no
  filter can undo — it turned out not to matter, because no released default-parser build carried
  the broken behaviour, but the shape of the risk is real. A future change that alters a value's
  content rather than a collection's membership would be unreproducible.

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
