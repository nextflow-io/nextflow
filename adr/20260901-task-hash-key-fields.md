# Recording the task hash key fields to explain cache misses from lineage

- Authors: Jorge Ejarque
- Status: draft
- Deciders: Jorge Ejarque, Paolo Di Tommaso
- Date: 2026-09-02
- Tags: lineage, task-hash, cache, resume, platform, nf-lineage

## Summary

Seqera Platform wants to answer "why was this task not cached?" by diffing two lineage `TaskRun`
records. That only works if the consumer knows which fields of the record actually feed the task
hash, and that set changes over time. This ADR frames how that field list could be produced,
versioned, and made available to the consumer, and records the options and their interactions. **No
option is selected yet** — see [Solution or decision outcome](#solution-or-decision-outcome).

## Problem Statement

`TaskHasher.compute()` builds the task hash from an **ordered, anonymous** `List<Object>` of key
values. The lineage `TaskRun` record (`nextflow.lineage.model.v1beta1.TaskRun`) records *most* of
the same information, so comparing the two records of two task runs is a plausible way to explain a
cache miss without re-running anything.

Three problems block that:

1. **The correspondence is implicit and currently wrong.** The record and the hasher were written
   and evolved independently, and they have already drifted (see
   [Current drift](#current-drift-between-taskhasher-and-the-lineage-record)). A consumer diffing
   the records today can report "no differences" for two tasks whose hashes genuinely differ.
2. **The correspondence is not one-to-one.** Some record fields are not hash inputs
   (`script`, `moduleId`, `workflowRun`), and `name` holds the *task* name (`FASTQC (3)`) while the
   hash uses the fully-qualified *process* name (`NFCORE_RNASEQ:RNASEQ:FASTQC`). A naive
   field-by-field diff produces false positives.
3. **The correspondence changes over time.** Keys are added to (or removed from) the hash as
   Nextflow evolves. A consumer reading a record written two years ago needs the field list that
   was true *then*, not the current one.

### Terminology

- **Task Hash Version (THV)** — an identifier for everything that determines the task hash bytes
  given the same inputs: the **set and order of the keys** `TaskHasher` feeds into the hash, the
  **per-key encoding** `HashBuilder` applies to each key value, and the **hash function** itself
  (`HashBuilder.DEFAULT_HASHING`, today `Hashing.murmur3_128()`). It changes only when one of those
  changes, so a single THV is shared by many Nextflow releases. It is **not** the lineage model
  version (`lineage/v1beta1`) and **not** the Nextflow version.

  The definition deliberately reaches beyond `TaskHasher`. Scoping the THV to the key list alone
  would leave it unchanged when `HashBuilder` changes how it encodes a `Path`, a collection, or a
  `Map` entry, or when the hash function is replaced — each of which moves every task hash while the
  field list stays identical. Under that narrower definition the comparison would report "no
  differences" for two tasks that genuinely differ, which is the exact failure this ADR exists to
  prevent. The consequence is that the THV cannot be owned by `TaskHasher` alone; it spans
  `TaskHasher` (keys) and `HashBuilder` (encoding, function), and those live in different modules
  (`nextflow` and `nf-commons`).
- **Hash-relevant field list** — for a given THV, the list of lineage `TaskRun` field names that
  contribute to the task hash. This is the artifact this ADR is about.

### There is no such identifier today

Nothing in the codebase currently identifies the task hash key set. A search for `hashVersion`,
`HASH_VERSION`, `CACHE_VERSION`, `SCHEMA_VERSION` and `hashSpec` across `modules/` and `plugins/`
returns nothing; `nextflow.cache` carries no version constant; and the term appears nowhere in
`adr/`, `specs/` or `docs/`. The THV is a new concept, and the closest existing constructs each fail
to serve as one:

| Existing construct | Location | Why it is not a THV |
| --- | --- | --- |
| `LinModel.VERSION` (`lineage/v1beta1`) | `nf-lineage/.../v1beta1/LinModel.groovy` | Versions the record *model*, not the hash. It is not bumped when `TaskHasher` gains a key while `TaskRun.groovy` is untouched — precisely the drift below. Also enforced by strict equality on read (`LinTypeAdapterFactory`), making it a compatibility gate rather than a history. |
| `BuildInfo.version` / `commitId` / `build` | `nf-commons/.../BuildInfo.groovy` | Too fine-grained: changes every release, whereas the THV changes rarely. Usable only as an inference source. |
| `WorkflowRun.metadata.nextflow` (`NextflowMeta.toJsonMap()`: `version`, `build`, `timestamp`, `enable`) | written in `LinObserver.storeWorkflowRun()` | Records the Nextflow version and the `enable` feature flags for every run. A useful secondary signal, but a property of the build, not a statement about the key set; inferring the key set from it is rejected below. |
| `CacheHelper.HashMode` | `nf-commons/.../CacheHelper.java` | A per-run runtime choice, orthogonal to the THV: two runs can share a THV and hash differently under `DEEP` or `SHA256`. Must be recorded separately (see prerequisites). |
| `HashBuilder.DEFAULT_HASHING` (`murmur3_128`) | `nf-commons/.../HashBuilder.java` | Not a version, but part of what the THV must cover — replacing it moves every hash with an unchanged key list. |

#### A candidate identifier exists in an open PR, but cannot be assumed

PR #6927 (`task-hasher-strategy`) would supply something very close to a THV: `TaskHasher` split into
an interface plus `TaskHasherV1`/`V2`/`V3` behind a `TaskHasherFactory.Version` enum whose values are
the strings `std/v1`, `std/v2`, `std/v3`, resolved once per session from `NXF_TASK_HASH_VER` and
cached on `Session.hashStrategy`. **It is an old PR — opened March 2026, last touched August 2026 —
and there is no guarantee it lands, so nothing below depends on it.** It is recorded here because it
would change the shape of several options, and because its existence is evidence about the problem
independently of whether it merges:

- Its versions **already disagree about the key set**: V3 hashes project `bin/` entries, V1 and V2 do
  not. So `binEntries` would be a hash-relevant field under one version and not another. A per-THV
  field list is not a hypothetical need.
- It would place the **per-key encoding** under the same version, via two new `HashBuilder` flags
  (`withOrderIndependentMaps`, `withCacheFunnelFirst`) configured by the strategy — narrowing, though
  not closing, the cross-module ownership question below. The hash function itself would stay an
  unversioned `nf-commons` constant.
- It would make the THV a **runtime choice rather than a build property**: `NXF_TASK_HASH_VER` lets
  two runs of the same Nextflow build hash under different versions. Any approach that infers the key
  set from the Nextflow version stops being sound at that point, which is one of the reasons it is
  not among the options.
- It would **not** supply the field list: keys are still collected into an anonymous `List<Object>`
  (`AbstractTaskHasher.collectKeys()`), so axis 1 remains open either way. A per-version
  `collectKeys()` override is, however, a natural shape for a per-THV list.

If it does not land, the THV has to be introduced by this work instead; if it lands later, the two
must be reconciled rather than coexisting as separate identifiers.

### Current drift between `TaskHasher` and the lineage record

| Hash key (`TaskHasher.compute()`) | Lineage `TaskRun` field |
| --- | --- |
| `session.uniqueId` | `sessionId` |
| `task.processor.name` (fully-qualified process name) | ⚠️ `name` holds the *task* name, a different string |
| `task.source` | ⚠️ only as `codeChecksum` (and `stubSource` under `-stub-run`) |
| container fingerprint | `container` |
| input name/value pairs | `input` |
| eval output commands | `eval` |
| task global vars (incl. `task.ext.*`) | `globalVars` |
| project `bin/` entries | `binEntries` |
| module resources bundle fingerprint | ❌ **absent** |
| `module` directive (environment modules) | ❌ **absent** |
| conda env / spack env / architecture | `conda`, `spack`, `architecture` |
| `'stub-run'` marker | ❌ **absent** |
| `HashMode` (`STANDARD`/`DEEP`/`LENIENT`/`SHA256`) | ⚠️ a `mode` *is* recorded, on every `Checksum` — but the wrong one (see below) |
| — | ➕ `script`, `moduleId`, `workflowRun` are recorded but are *not* hash inputs |

The three absent entries are the failure mode that matters: with identical field values, an upgrade
that changed the `module` directive or a `-stub-run` execution produces a different hash while the
record diff shows nothing.

`HashMode` is a near-miss rather than an absence, and the distinction is worth stating precisely
because it changes what remains to be done. `Checksum` already carries a `mode` field, and
`TaskRun.codeChecksum`, `binEntries` and the input paths all use it, so a consumer does see a mode.
It is not the one that matters:

- **Wrong scope.** `Checksum.ofNextflow()` stamps `CacheHelper.HashMode.DEFAULT()` — the process-wide
  default, set once from `NXF_CACHE_MODE`. `TaskHasher.compute()` hashes under
  `task.processor.getConfig().getHashMode()`, which is `HashMode.of(config.cache) ?: DEFAULT()` — a
  **per-process** value driven by the `cache` directive. A pipeline where one process declares
  `cache 'deep'` and the rest do not records a single identical mode for every task, and it is the
  wrong one for exactly the process whose caching behaviour differs.
- **Wrong subject.** The field describes how *that checksum* was computed, not how the task hash was.
  Reading it as the latter is an inference the model does not license, and one that silently becomes
  false as soon as either side changes.
- **Historically mislabelled.** Until #7582 (`2c68fa1c9`), `Checksum.ofNextflow()` computed its value
  via `CacheHelper.hasher(value)` — hard-wired to `HashMode.STANDARD` — while labelling the result
  with `DEFAULT()`, so under `NXF_CACHE_MODE=DEEP` the stored checksum was a standard one carrying
  the label `deep`. That is fixed on master, but records written by earlier releases carry the wrong
  label and cannot be trusted retroactively even for their own stated purpose.

So the prerequisite below narrows but does not disappear: what is missing is not "a mode" but *the
effective per-process mode `TaskHasher` used*, recorded as a property of the task run rather than of
one of its checksums.

## Goals or Decision Drivers

- **Honesty over completeness.** The comparison is *explanatory*, not *verifiable* (see Non-goals).
  It may say "I cannot see this key", but it must never assert "no differences" when the hashes
  differ.
- **Resistance to drift.** The mechanism must make it structurally hard — ideally impossible — for
  the field list to fall out of sync with `TaskHasher`. The drift table above is evidence that
  convention plus code review does not achieve this.
- **Resolvable for records already written.** Every lineage record in existence today carries no
  THV. The design must still be able to interpret them.
- **Resolvable for records from unknown builds.** Platform will encounter edge releases, dev builds,
  and versions newer than itself. An unresolvable THV must be rare or impossible, not routine.
- **No cache invalidation.** Nothing in this work may change the bytes fed to the hasher; an
  accidental change silently invalidates every user's cache on upgrade.
- **Negligible cost per task.** Lineage stores are dominated by `TaskRun` objects; per-task overhead
  multiplies by the task count, which can reach 10^5–10^6.

## Non-goals

- **Verifiable attribution.** This ADR deliberately targets an *explanatory* diff ("these fields
  differ"), not a proof that the reported set of differing keys is exactly complete. Proving
  completeness would require recording every key's individual contribution plus the hash mode and
  re-deriving the hash, which is a strictly larger design and a separate decision.
- **Recomputing task hashes** outside Nextflow, in Platform or anywhere else.
- **The Platform-side data model and UI** for presenting the diff.
- **Changing the task hash itself**, its key set, or its ordering.
- **Explaining non-hash cache misses** — a missing/corrupt work directory, a failed exit status, or
  a cache DB miss are separate causes and out of scope here.

## Considered Options

### Two payloads, not one

The options below move two *different* pieces of information. Treating them as one makes the axes
look independent and leaves "what exactly gets embedded?" unanswered, so they are separated here
before the options are listed.

| | Size | Answers | Changes when |
| --- | --- | --- | --- |
| **Field list** | ~200 bytes | *"Which record fields do I diff, and which differences are irrelevant?"* | the key set changes |
| **THV** | a few bytes | *"Are these two records comparable at all, and is a no-difference diff trustworthy?"* | the key set, the per-key encoding, **or** the hash function changes |

They are not two encodings of the same fact, and the relation between them runs one way only: a
different field list implies a different THV, but **an identical field list does not imply an
identical THV**. The list names *lineage record fields* — it says *what* was hashed. The per-key
encoding and the hash function decide *how*, and no list of field names can express either.

The clearest illustration is in PR #6927, which is not assumed to land but does demonstrate that the
case is real rather than theoretical. Its `TaskHasherV1` and `TaskHasherV2` both inherit
`collectKeys()` unchanged from `AbstractTaskHasher` and differ *only* in `createHashBuilder()`:

| | key set | `withOrderIndependentMaps` | `withCacheFunnelFirst` |
| --- | --- | --- | --- |
| `TaskHasherV1` | inherited | `false` | `false` |
| `TaskHasherV2` | inherited — **identical to V1** | `true` | `true` |
| `TaskHasherV3` | overrides `collectKeys()` — adds `bin/` entries | `true` | `true` |

V2 → V3 is the case intuition catches: the key set changes, so the field list changes with it. V1 → V2
is the case it misses: identical field list, identical record values, identical task sources,
**different hash**, because a `Map` is serialised by its values under V1 and by its `entrySet()` under
V2. A consumer holding the list but not the THV compares those two records, finds nothing different,
and reports "no differences" for two tasks whose hashes genuinely differ — the precise failure the
honesty driver forbids, reached without anyone making a mistake.

The THV alone is not sufficient either: it says *whether* to trust a diff, not *what* to diff.

**Consequence: the THV is recorded under every option, without exception** — including the options
that embed the list, where it is not redundant but is the only witness to the encoding and the hash
function. Only the *field list's* location is genuinely open. It also means the axis-2 entries are
not mutually exclusive alternatives to be picked between; they are positions on a distance ladder,
and a workable answer may well combine two of them.

**Axis 1 — how the field list and the THV stay true to `TaskHasher`:**

- A1. Convention (a comment saying "bump the version if you change these keys")
- A2. A guard test asserting the ordered key names for a canonical task
- A3. Derive the list from `TaskHasher` itself, so it cannot be stated separately

**Axis 2 — how far the field list sits from the record it describes:**

| | Option | What it stores | Distance from a `TaskRun` | THV still required? |
| --- | --- | --- | --- | --- |
| B1 | In every `TaskRun` | list **and** THV, per task | none | yes — the list covers the key set, the THV covers encoding + hash function |
| B2 | Once per run, in `WorkflowRun` | list + THV per run, THV also stamped per task | one hop, already traversed | yes, same reason; the extra `TaskRun` stamp is an early exit, not a capability |
| B3 | Sibling metadata file, one per store | list, keyed by THV | store layout; may be unreachable | yes — stamped on the record; **not standalone** |
| B4 | Central registry, one for all stores | list, keyed by THV, outside the records | out of band | yes — stamped on the record; **not standalone**, the registry has no other entry point |

Read that way, **B3 and B4 hold the same artifact** — a THV → list map — and differ in who owns it and
how far it reaches: B3 keeps a per-store copy written by the run, B4 keeps one global copy written by
a release. **Neither is implementable alone**: a map needs a key, so both require the THV stamped on
the record, which is B1 or B2 in miniature.

A fifth possibility — store nothing and have the consumer infer the key set from the Nextflow version
already in `WorkflowRun.metadata.nextflow`, via a version → THV table held outside Nextflow — is
deliberately **not** listed as an option. It is a guess about a build rather than a statement by it,
so a patched or vendored Nextflow reports a version whose key set may not be the one it used; dev,
edge and snapshot builds do not map onto version ranges, and those are disproportionately the builds
whose users are debugging a cache miss after an upgrade; and a version-keyed lookup always returns
*some* row, so it fails confidently rather than admitting a gap, which the honesty driver rules out.
Records written before any stamp exists are covered instead by B4's `null` → THV_0 entry.

## Pros and Cons of the Options

### Axis 1

#### A1 — Convention

*What it means.* A comment above the key accumulation in `TaskHasher.compute()` saying "if you add
or remove a key here, update the field list and bump the THV", and a hand-written field list living
wherever it is stored. Nothing enforces either step; correctness rests on the author noticing and on
review catching it if they do not. This is what the codebase does today for the hasher/record
correspondence, minus the comment.

- Good, because it costs nothing to implement.
- Bad, because it is the status quo, and the status quo has already failed three times
  (module bundle fingerprint, `module` directive, `stub-run`).
- Bad, because the failure is silent and lands in a user-facing answer. A stale list does not
  degrade the feature, it makes it lie.

#### A2 — Guard test

*What it means.* A Spock test builds a canonical `TaskRun` fixture, asks `TaskHasher` for the ordered
list of key names it would feed the hash, and asserts it equals a literal list checked into the test.
Adding, removing or reordering a key makes the build red until the author edits that literal — at
which point the field list and the THV constant, kept adjacent to it, are in front of them. It
enforces *attention at the right moment*, not correctness of what they then write.

- Good, because it makes any change to the hash key set fail the build, forcing the author to touch
  the field list in the same commit.
- Good, because it is cheap: one Spock test with a fixture listing the expected ordered key names.
- Good, because it doubles as a regression test against accidental hash changes — valuable on its
  own merits, independent of this feature.
- Bad, because the author can still update the fixture without updating the field list; it enforces
  *attention*, not *correctness*.

#### A3 — Derive from `TaskHasher`

*What it means.* Refactor the anonymous `keys << value` accumulation into a named-key builder
(`keys.put('condaEnv', conda)`), so the ordered name list becomes a property of the code and the
field list is generated rather than written. The *key-set component* of the THV then follows from
that name list by content derivation; the encoding and hash-function components do not, and still
have to be supplied from the `HashBuilder` side.

- Good, because drift becomes structurally impossible: changing a key changes the generated list and
  the THV, with no human step to forget.
- Good, because it establishes a **single canonical vocabulary** shared by `TaskHasher`, the lineage
  `TaskRun` model, and `-dump-hashes`, which collapses problem 2 (the record-field mapping) into
  naming rather than a maintained table.
- Good, because it fixes `-dump-hashes` as a side effect: `dumpHashEntriesJson()` currently emits
  only `hash`/`type`/`value` per entry with no name, so today's primary debugging tool for cache
  misses produces an unlabelled list. Named keys make it directly readable.
- Bad, because it is a refactor of a cache-critical class. The key *order* and the per-key encoding
  must be proven byte-identical before and after, or every user's cache is invalidated on upgrade.
- Bad, because conditional keys (`container` only when enabled, `spack`/`arch` only when set) mean
  the emitted list is per-task, not static — the name list must be the *declared* set, not the
  *emitted* set, or the THV becomes task-dependent.

### Axis 2

#### B1 — Embedded in every `TaskRun` record

*What it means.* Every `TaskRun` record carries two new fields: the THV, and the full ordered field
list for it. There is no THV → list map anywhere, because each record answers for itself; the map
exists only implicitly, replicated once per task. A consumer needs nothing but the two records it is
already diffing.

- Good, because every record is fully self-describing; no resolution step, no lookup, no unknown-THV
  case, ever.
- Bad, because the cost multiplies by task count: ~14 names is ≈200 bytes of JSON, so ≈20 MB of pure
  duplication on a 100k-task run, in the store's most numerous object type.
- Bad, because it is derived data replicated 100k times, which invites divergence between copies if
  anything ever writes records through a different path.
- Bad, because — as with B2 — it can only describe records written after the fields exist. Every
  record in the store today predates them, so a legacy path is needed either way; this is the real
  limitation of both embedding options, and the only one the external tables answer natively.

#### B2 — Embedded once per run, in the `WorkflowRun` record

*What it means.* The `WorkflowRun` record carries the THV and its field list; each `TaskRun` carries
the THV only. The THV → list map is degenerate — exactly one entry, stored inside the run it
describes. A consumer reads `TaskRun.taskHashVersion` to decide whether two tasks are comparable, and
follows the existing `TaskRun.workflowRun` link when it needs the list itself.

- Good, because it keeps all of B1's self-description with none of its size cost — one copy per run.
- Good, because the resolution path already exists: `TaskRun.workflowRun` points at the
  `WorkflowRun` record, which is where the Nextflow version already lives
  (`WorkflowRun.metadata.nextflow`, populated in `LinObserver.storeWorkflowRun()`). One extra fetch
  that consumers already perform.
- Good, because it works offline, air-gapped, and for dev/edge builds Platform has never heard of.
- Good, because it requires no cross-service release coordination: a Nextflow release that changes
  the key set ships its own explanation.
- Bad, because a wrong list is **immutable** — if a release ships a bad derivation, every record it
  wrote is wrong forever and the diff misreports with confidence.
- Bad, because it does not help records already written (no such field exists today).

#### B3 — Sibling metadata file in the lineage store / cache directory

*What it means.* Nextflow writes a THV → list map as a plain file in the lineage store (say
`.meta/task-hash-versions.json`, beside the record tree), adding an entry the first time it writes a
record under a THV that file does not yet have. The map is therefore **local to one store** and
covers only the THVs that store happens to contain. Records still carry the THV stamp — without it
nothing selects an entry from the file — so B3 is *B1/B2's stamp plus an out-of-record list*, not an
alternative to them. Its only genuine saving over B2 is that the `WorkflowRun` model is untouched.

- Good, because it leaves the `WorkflowRun` model untouched and adds only the THV stamp to
  `TaskRun` — the smallest record-model change of any option that keeps the list reachable offline.
- Bad, because reachability is unproven: if Platform ingests lineage through the API or the
  `nf-tower` plugin rather than listing the store, a loose file in the store is invisible to it.
- Bad, because it invents a second, unversioned, unmodelled artifact in a store whose whole point is
  that everything in it is a typed, versioned record.
- Bad, because it does not survive record export, copy, or partial ingestion — the record and its
  explanation can be separated.
- Bad, because the map is **mutable shared state**: every run writing into the store may have to
  update it, and updating means read-modify-write. Two concurrent runs each add their entry and the
  second write drops the first, leaving records whose stamped THV resolves to nothing. On
  object-store-backed lineage stores (S3, GCS) last-writer-wins is the *only* available semantics —
  no atomic update, no append, no lock — so this is the normal case rather than a narrow race, and
  it gets likelier exactly when a site is running mixed Nextflow versions against a shared store,
  which is the situation the map exists to describe.
- Bad, because loss is not the worst outcome. One slot per THV means the last writer also decides the
  *content* of an existing entry, so one buggy, patched or vendored build writing a different list
  under the same THV silently reinterprets every record already written under it — the
  "confidently wrong" failure, now reachable without anyone shipping a bad release.
- *Mitigation, for completeness.* One write-once file per THV (`.meta/thv/<id>.json`, never
  rewritten) removes the concurrency hazard — but by then B3 has grown into a content-addressed
  sidecar store of unmodelled files, markedly more machinery than the single `WorkflowRun` field
  B2 costs.

#### B4 — Central registry (THV → field list)

*What it means.* One THV → list document, maintained outside every store and covering **every THV
that has ever existed**, not just the ones a given store contains. It is versioned and released like
code — shipped as a resource in the Nextflow distribution, mirrored or hosted by Platform, or both —
and is edited by a release process rather than written at run time. Records still carry the THV
stamp, so B4 likewise is not standalone.

*Difference from B3.* Same data structure, different owner, scope and lifecycle: B3's copy is written
by the run, travels with the data, and knows only what that store saw; B4's copy is written by a
release, is reachable independently of any store, and knows the whole history. The practical
consequences follow from that — B3 can never explain a THV its store did not produce, and B4 can be
corrected once for every record ever written.

- Good, because it is **correctable**: a wrong list can be patched centrally and every historical
  record is retroactively fixed. This is the one thing embedding cannot do, and given the drift
  already found, "we shipped a wrong list once" is not hypothetical.
- Good, because it stores the fact once, no matter how many runs or records exist.
- Good, because it can carry richer per-THV information later (deprecation notes, per-key
  descriptions, human-readable explanations) without touching the record model.
- Good, because **an absent THV is itself a usable key**: the registry can hold a `null` entry
  pointing at a designated legacy THV_0, so every record written before the stamp existed resolves
  through the same lookup as any other, with no new record field and no version arithmetic. This is
  the only mechanism among the options that addresses records already written, and it is why B4
  cannot simply be dropped in favour of embedding. Note the limit — one default collapses every
  historical key set into a single answer, so it is honest only if THV_0's list is the conservative
  intersection of the key sets that window spans, or is explicitly marked "best effort, unverified".
- Bad, because it makes every lookup a dependency on an artifact Platform must have, for a THV it may
  never have seen — dev builds, edge releases, and any Nextflow newer than the deployed Platform all
  land in the unresolvable case. Under the honesty driver that degrades to "I don't know", which is
  precisely the outcome the feature exists to avoid.
- Bad, because it needs the mapping retained for every THV that ever wrote a record, forever, plus a
  release process step to publish new entries — a maintained-forever table for a fact the record
  could simply carry.
- Bad, because it splits the truth across two release trains (Nextflow and Platform) that ship at
  different cadences.

### How the axes interact

They are separable but not independent: the axis-2 choice determines how much correctness the axis-1
choice has to supply, because it determines *when* the list is produced and *who* can still fix it.

| | B1 / B2 — embedded | B3 — sibling file | B4 — central registry |
| --- | --- | --- | --- |
| **A1** convention | Worst cell. A wrong list is frozen into records immutably, and a missed bump is undetectable by anyone. | Same, and the file can be separated from the records it explains. | Wrong in a subtler way: records claim a THV whose content silently moved, so the registry answers *confidently wrong*. |
| **A2** guard test | Viable. The build fails on any key change, so the embedded list is at least *attended to* every time it could rot. | Viable, reachability aside. | Viable — A2 is precisely what makes a hand-bumped THV trustworthy enough to key a registry with. |
| **A3** derived | Natural fit. The list is free at write time and cannot be stale. | Pays B3's reachability cost for no gain. | Works, but the registry's correctability advantage buys much less once the list cannot be wrong. |

Four concrete dependencies behind that table:

1. **Embedding presumes the list can be produced at runtime.** B1/B2 write whatever the running
   Nextflow believes the list to be. Under A3 that belief is derived and correct by construction;
   under A1/A2 it is a human-written constant frozen into every record that release writes. B2's
   "a wrong list is immutable" con is therefore largely an *A1/A2* con — A3 mostly retires it, and
   an embedding option is correspondingly harder to justify without A3 behind it.
2. **A registry presumes the THV is trustworthy.** B4 resolves the list *by identifier*, so a
   THV that failed to move when the keys moved returns a wrong list with full confidence and no
   detectable symptom — strictly worse than an unresolvable THV, which at least degrades to "I don't
   know". A derived (A3) or test-pinned (A2) THV is a precondition for B4, not a refinement of it.
3. **A3 changes the payload's shape, not just its provenance.** With a single vocabulary shared by
   `TaskHasher` and the lineage record, the field list is a list of names. Without it, what must be
   stored is a *translation table* — hasher key → record field, plus an explicit "unrepresented" set
   — roughly double the size and itself a hand-maintained artifact, which is the very thing this ADR
   exists to stop creating. The ~200 byte figure quoted in B1/B2 assumes A3.
4. **Axis 1 cannot help records already written.** Whatever mechanism is adopted, it starts working
   on the day it ships; every record in the store today predates it. Only B4's `null` → THV_0 entry
   addresses those, and populating it is archaeology of `TaskHasher`'s git history (and of
   `HashBuilder`'s, for the encoding component) rather than anything a forward-looking mechanism can
   derive. That work is the same size whichever axis-1 option is chosen, and it is worth scoping
   separately from the rest.

## Solution or decision outcome

**Not decided yet.** This ADR is a draft: the options and their interactions are recorded above so
the choice can be made deliberately, not reached by default.

What still has to be settled, in rough dependency order:

1. **Axis 1** — whether the field list is derived from `TaskHasher` (A3), pinned by a guard test over
   a hand-written list (A2), left to convention (A1), or some combination. A3 is a refactor of a
   cache-critical class, so it cannot be committed to until a test has demonstrated that the hash
   produced for a fixture task is byte-identical before and after the rename.
2. **Axis 2** — where the field list lives, bearing in mind that the entries are positions on a
   ladder rather than exclusive alternatives, that B3 and B4 both need a stamp to be usable at all,
   and that only B4 addresses records already written.
3. **Whether the THV is introduced by this work or inherited.** If PR #6927 lands it supplies one;
   if it does not, this work has to define it. The two must not both happen.
4. **The record-completeness prerequisites** below, which are needed under every combination and are
   the part that can start independently of the rest.

Whatever is chosen has to satisfy the drivers above — in particular it must never report "no
differences" for two tasks whose hashes differ, and it must produce an answer when the two runs were
hashed under different THVs, since that is the common case rather than an edge one.

## Prerequisites, whichever options are chosen

None of the options above fixes the drift; they only describe it. Under the honesty driver, the
following must land alongside this work, or the diff will report "no differences" on tasks with
different hashes:

1. Add the missing hash inputs to the lineage `TaskRun` record: the `module` directive values, the
   module resources bundle fingerprint, and the `stub-run` marker. For `HashMode`, record the
   *effective per-process* mode (`task.processor.getConfig().getHashMode()`) as a field of the task
   run — the existing `Checksum.mode` is the global `NXF_CACHE_MODE` default and describes the
   checksum it sits on, not the task hash, so it cannot stand in for it.
2. Record the fully-qualified process name used by the hasher, distinctly from the tagged task
   `name` — or exclude `name` from the field list, since it is not a hash input.
3. Ensure a new hash key cannot be added silently: every key the hasher feeds must be either present
   in the field list or explicitly declared unrepresented. Under A2 or A3 the guard test is the place
   to assert this; under A1 there is no such place, which is itself an argument on axis 1.

`HashMode` deserves emphasis: it is a *runtime* choice, not part of the THV. Two runs can share a THV
and an identical field list and still hash differently because one ran under `DEEP` or `SHA256`, and
any mode added later widens that gap rather than closing it. The `mode` already present on
every `Checksum` does not cover this — it is the global default rather than the per-process effective
mode, and it describes a checksum rather than the task hash — and in records written before #7582 it
is mislabelled as well. The cost of doing it properly is one enum-valued field per `TaskRun`, which
is the cheapest item on this list and guards a whole class of otherwise invisible miss.

## Open questions

- **Lineage model version.** Adding fields to `TaskRun`/`WorkflowRun` is read-compatible in both
  directions with Gson (new readers see `null`, old readers ignore unknown fields), and
  `LinTypeAdapterFactory` gates on strict equality of `version` (`lineage/v1beta1`). So no bump is
  strictly required — but whether growing the model without a bump is acceptable practice should be
  settled, since `WorkflowRun` also feeds the execution hash
  (`CacheHelper.hasher(value).hash()` in `LinObserver.storeWorkflowRun()`) and new fields will change
  that hash.
- **THV form.** A monotonic integer is friendlier in a UI; a content-derived short hash cannot be
  forgotten or collide. A hybrid (integer, with the guard test asserting it was bumped when the
  derived content hash changed) gets both. If PR #6927 lands, its `std/vN` enum value is the obvious
  candidate and this question is answered for us — but it is a *name*, and names can be rebound: an
  edited `std/v1` would stamp a version string that no longer means what it did, so it would still
  need fixture tests pinning each strategy.
- **THV ownership across modules.** Because the THV covers encoding and the hash function as well as
  the key list, it cannot be derived from `TaskHasher` alone — `HashBuilder` lives in `nf-commons`,
  which has no visibility of `TaskHasher`. Options: derive the key-list part in `nextflow` and fold
  in an encoding/function component declared in `nf-commons`; or keep a single hand-bumped THV
  constant in `nf-commons` guarded by fixture tests on both sides. PR #6927 would narrow this by
  having each strategy configure `HashBuilder` itself, leaving only the hash function
  (`DEFAULT_HASHING`) unversioned — but only if it lands, and the residue would still need pinning.
- **Whether to wait on PR #6927 at all.** It has been open since March 2026. Building on it risks
  blocking indefinitely; building beside it risks two competing hash-version identifiers that then
  have to be merged. Worth resolving with its author before any implementation starts, since it
  determines whether the THV is this ADR's to define.
- **Whether Platform can read the lineage store directly.** It determines whether B3 is feasible at
  all; worth confirming so the option can be closed on evidence rather than left open.
- **`-dump-hashes` alignment.** Whether the named-key output should be promoted to a documented,
  stable format consumable by tooling, or remain a debug aid.

## Links

- Refines [Data lineage](20250508-data-lineage.md)
- Related PR #6927 "Add versioned task hasher strategy" (`task-hasher-strategy`, open since March 2026,
  **not assumed to land**) — would introduce `TaskHasherFactory.Version` (`std/v1..v3`) and
  `NXF_TASK_HASH_VER`, and carries the module resources bundle fingerprint into the hash (#6914, closes
  #6128)
- Related PR #7171 "ADR: lineage record version compatibility" (open) — covers the lineage model
  version question raised in the open questions above
- Fixed upstream: #7582 (`2c68fa1c9`) "Fix lineage checksum computed in standard mode regardless of
  label" — `Checksum.ofNextflow()` now hashes under the mode it labels
- Code: `modules/nextflow/src/main/groovy/nextflow/processor/TaskHasher.groovy`,
  `modules/nf-lineage/src/main/nextflow/lineage/model/v1beta1/TaskRun.groovy`,
  `modules/nf-lineage/src/main/nextflow/lineage/LinObserver.groovy`
