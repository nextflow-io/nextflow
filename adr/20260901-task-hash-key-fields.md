# Recording the task hash key fields to explain cache misses from lineage

- Authors: Jorge Ejarque
- Status: proposed
- Deciders: Ben Sherman, Paolo Di Tommaso
- Date: 2026-09-02
- Updated: 2026-10-05
- Tags: lineage, task-hash, cache, resume, platform, nf-lineage

## Summary

Seqera Platform wants to answer "why was this task not cached?" by diffing two lineage `TaskRun`
records. That only works if the consumer knows which fields of the record actually feed the task
hash, and that set changes over time. This ADR frames how that field list could be produced,
versioned, and made available to the consumer, and records the options and their interactions.

**The proposal is A2 + B4**: each published hash version is a declarative spec file pinned by guard
tests, every `TaskRun` record carries the id of the version that hashed it, and the specs are served
from a central registry that also maps each key to the record fields a consumer should diff. See
[Solution or decision outcome](#solution-or-decision-outcome).

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

PR #6927 (`task-hasher-strategy`) would supply something close to a THV: `TaskHasher` split behind a
`TaskHasherFactory.Version` enum whose values are `std/v1`, `std/v2`, `std/v3`. **It has been open
since March 2026 and may never land, so nothing below depends on it.** It is recorded because it is
evidence about the problem regardless of whether it merges:

- Its versions **already disagree about the key set** — V3 hashes project `bin/` entries, V1 and V2
  do not — so a per-THV field list is not a hypothetical need.
- Its V1 and V2 have an **identical key set and still hash differently**, because the strategy
  configures `HashBuilder` flags. A field list alone cannot express that, which is the basis of the
  "two payloads" argument below.
- It would **not** supply the field list: keys are still collected into an anonymous `List<Object>`,
  so axis 1 stays open either way.

The `TaskHasherFactory` extension point has since landed on master on its own, in #7638
(`631ed9734`), and carries no version identifier. So a versioned hasher plugs in without #6927, and
the decision below introduces the THV rather than inheriting it.

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

It returns in the decision below, but only in a narrowed form and only where no stamp can ever exist:
as an explicitly *inferred* answer for records already written, drawn from a table of exact release
tags rather than open version ranges, returning no answer at all for a release not in the table. That
is a different claim from the one rejected here. The rejection stands for every record a stamped
Nextflow writes, which is every record from now on.

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
- Good, because **an absent THV can still be resolved**: the registry can serve a second small
  artifact mapping Nextflow releases to THVs, so a record written before the stamp existed is
  answered without any new record field. This is the only mechanism among the options that addresses
  records already written, and it is why B4 cannot simply be dropped in favour of embedding. Note the
  limit — the answer is an inference about the build rather than a statement by it, so it has to be
  presented as one, and releases older than the first reconstructed THV must resolve to nothing. The
  decision below takes this shape; an earlier draft used a single `null` → legacy entry instead,
  which collapsed all of history into one answer and was dropped.
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
   on the day it ships; every record in the store today predates it. Only B4 addresses those, through
   a release → THV table served beside the specs, and populating it is archaeology of `TaskHasher`'s
   git history (and of `HashBuilder`'s, for the encoding component) rather than anything a
   forward-looking mechanism can derive. That work is the same size whichever axis-1 option is chosen, and it is worth scoping
   separately from the rest.

## Solution or decision outcome

**Proposed: A2 + B4.**

- **A2** — each hash version is a declarative spec file in the Nextflow jar, pinned by guard tests
  that compare real hashes rather than key names. A3 was not chosen: deriving the list from
  `TaskHasher` cannot express the per-key encoding or the hash function, which the "two payloads"
  analysis above shows is the half that fails silently.
- **B4** — the specs are served from a central registry, keyed by an identifier stamped on every
  `TaskRun` record. The field mapping travels inside the spec file, so a consumer holding a stamp
  fetches exactly one document.

### 1. A hash version is a spec file

A **task hash version** (THV) is a `TaskHashSpec`: an ordered list of keys with the contributor that
produces each key's value, a set of encoding rules, and a hash function. It is loaded from a JSON
resource rather than written in Groovy, so a version is data that can be published, diffed and served
unchanged, not code that has to be kept compiling.

Identifiers are `std/<major>.<minor>`. The minor number moves when the current code can still
reproduce the older version; the major number moves only when it cannot, which is the honest signal
that an old cache has become unexplainable. Seven versions are published, covering every release from
25.10.0 to today:

| THV | First release | Ended by | What changed |
| --- | --- | --- | --- |
| `std/v1.1` | 25.10.0 | #6605 `1ca327c80` | `isAssetFile` widened to match the asset root, not only `baseDir` |
| `std/v1.2` | 25.11.0-edge | #6643 `295f17307` | the v2 parser became the default and stopped collecting `params.*` refs |
| `std/v1.3` | 26.01.1-edge | #6679 `d54ff29af` | record-types encoding: order-independent maps, funnel first |
| `std/v1.4` | 26.03.0-edge | #7165 `785e801ad` | `params.*` refs folded back into the task global vars |
| `std/v1.5` | 26.04.2 | #6914 `029e52eef` | the module resources bundle became a key |
| `std/v1.6` | 26.08.0-edge | #7575 `9e7a492a3` | eval outputs hashed as the raw map, not the derived string |
| `std/v1.7` | 26.09.0-edge | — | current master |

Two findings from building this are worth recording, because both contradict what the option analysis
assumed:

1. **Three encoding rules, not two.** PR #6927 carries `withOrderIndependentMaps` and
   `withCacheFunnelFirst`. Reproducing 25.10.0 needs a third, `withAssetRootDetection`, because #6605
   changed which files are hashed by content rather than by metadata. Nothing about that change is
   visible in the key set, which is the "two payloads" argument meeting a real case.
2. **The encoding rules leaked.** `HashBuilder` applied them to top-level values only. A value nested
   in a collection, or reached through a `CacheFunnel`, was hashed by a fresh builder carrying the
   defaults, so a historical rule set stopped applying exactly where it mattered. `std/v1.1` and
   `std/v1.2` produced identical hashes until this was fixed. The same leak is present in #6927.

### 2. Keeping the versions honest (A2)

`StdSpecsGoldenTest` holds the whole of axis 1. It builds one task that fires all thirteen keys and
asserts three things:

1. **The newest spec reproduces the default hashing path.** This is the drift guard the ADR asked
   for. The moment someone changes what master hashes, the newest spec stops describing it and the
   build fails, which is precisely when a new spec is due. The test also asserts that no key is
   silently empty, otherwise the comparison would pass vacuously.
2. **Every published spec still produces the hash it was published with.** The expected values are
   frozen. A moved value means an already published version changed, so every cache written under it
   has stopped resolving.
3. **Every published spec still has its published fingerprint** — a hash of its id, ordered keys,
   contributor names, encoding rules and hash function. This catches a changed definition that the
   task hash cannot see: `std/v1.1` and `std/v1.2` differ only in asset detection, which needs a file
   inside a Git repository and so hashes identically in a unit fixture.

Guard tests alone would only prove the specs are self-consistent. Each version was therefore also run
against the genuine release it claims to describe: a pipeline exercising all thirteen keys was run on
that release, then resumed under the matching THV. All seven resumed 8 of 8 tasks
`CACHED`, and the default path did the same against 26.09.0-edge. A cross-check confirmed the specs
are mutually distinct in the predicted places rather than accidentally equal.

### 3. The stamp

`TaskRun.hashVersion` records the id of the version that produced `TaskRun.hash`. `AgentRun` carries
the same field. A hasher reports its own version; the default path reports `StdSpecs.latest().id`,
which is sound precisely because guard 1 above fails the build if the two ever diverge.

This is the B1-in-miniature that B4 requires. It is one short string per task, and it is the only
witness to the encoding rules and the hash function, neither of which any list of field names can
express.

### 4. The field mapping lives in the spec file

Each key in a spec gains a `lineage` list naming the `TaskRun` fields a consumer should diff for that
key. An empty list means the key is not recorded, and the consumer must say so rather than report no
difference.

```json
{ "key": "SPACK", "contributor": "spackEnvAndArch", "lineage": ["spack", "architecture"] }
```

Putting the mapping in the spec rather than in a second document means a consumer holding a stamp
makes exactly one request and gets everything that stamp implies — what was hashed, in what order,
under which encoding, and where to look for each key in the record. It also means the mapping is
versioned with the thing it describes, by construction.

The full mapping, as of `std/v1.7`:

| Hash key | Lineage field(s) | Note |
| --- | --- | --- |
| `SESSION_ID` | `sessionId` | |
| `PROCESS_NAME` | `name` | superset, see below |
| `TASK_SOURCE` | `codeChecksum` | |
| `CONTAINER` | `container` | |
| `INPUTS` | `input` | |
| `EVAL_OUTPUTS` | `eval` | |
| `SCRIPT_VARS` | `globalVars` | |
| `BIN_ENTRIES` | `binEntries` | |
| `RESOURCES_BUNDLE` | `resourcesBundle` | **new field** |
| `ENV_MODULES` | `envModules` | **new field** |
| `CONDA` | `conda` | |
| `SPACK` | `spack`, `architecture` | one key, two fields |
| `STUB_MARKER` | `codeChecksum` | shares a field with `TASK_SOURCE` |

Only two lineage fields are added: `resourcesBundle` (String, the bundle fingerprint) and
`envModules` (`List<String>`, the `module` directive values). The ADR's prerequisite list asked for a
third, a separate field for the stub-run marker; it is not needed, because a stub run already changes
`codeChecksum` — the stub block is the source that gets hashed. Pointing `STUB_MARKER` at
`codeChecksum` is therefore accurate and costs nothing.

`PROCESS_NAME` → `name` needs a caveat in the consumer. The hasher uses the fully-qualified process
name; `name` is the tagged task name, which contains it plus the tag. So `name` differing does not
prove `PROCESS_NAME` differed, but `name` being equal does prove it did not. The diff may over-report
this key and can never under-report it, which is the direction the honesty driver allows. The earlier
prerequisite to add a separate process-name field is dropped: a second field would be cheaper to read
but is not worth a per-task cost for a key that practically never differs between two records being
compared.

The key for the resources bundle is named `RESOURCES_BUNDLE`, not `MODULE_BUNDLE`, because three
unrelated things around it are already called "module": the remote Nextflow module in `moduleId`, the
environment modules in `ENV_MODULES`, and the `resources/` bundle itself. `ResourcesBundle` is also
the class name in the code. Key names are part of the published contract once a spec reaches the
registry, so they cannot be corrected afterwards without a new version.

### 5. The registry (B4)

Published under `nextflow-io/schemas`, beside the other Nextflow schemas and reachable the same way:

```
task-hash/v1/
  schema.json              # the shape of a spec file
  legacy-releases.json     # Nextflow release -> THV, for records with no stamp
  specs/std-v1.1.json      # byte-identical to the jar resource
  ...
  specs/std-v1.7.json
```

There is no index of available specs, because nothing needs one: a stamped record names its spec, so
the consumer fetches `specs/std-v1.7.json` by that name. A spec file is immutable once published and
can be cached forever.

Two guards keep the registry and the jar in step, and they are the part that makes B4 safe rather
than merely convenient:

1. **At release time**, publishing a new spec copies the jar resource byte for byte. The release
   fails if a spec file that already exists in the registry differs from the jar.
2. **On a schedule**, a CI job compares every published spec against the jar resource of the same id
   and fails on any difference. This catches an edit made directly in the registry, which is the one
   failure mode B4 adds over the embedding options: records claiming a THV whose content silently
   moved, answered confidently and wrongly.

### 6. Records written before the stamp

Every record in existence today has no `hashVersion`, and nothing can be added to them retroactively.
They are resolved from the Nextflow version already present in the run, reached via
`TaskRun.workflowRun` → `WorkflowRun.metadata.nextflow.version`, against `legacy-releases.json`:

```json
{
  "formatVersion": "task-hash/v1",
  "unidentifiableBefore": "25.10.0",
  "versions": [
    { "id": "std/v1.1", "nextflow": ["25.10.0", "25.10.1", "...", "25.10.8"] },
    { "id": "std/v1.7", "nextflow": ["26.09.0-edge", "26.09.1-edge"] }
  ]
}
```

Four properties make this acceptable where the rejected fifth option was not:

- **It applies only where no stamp can exist.** A record written by any Nextflow that carries the
  stamp is resolved from the stamp. This table never overrides one.
- **The answer is labelled inferred, not asserted.** The objection that a patched or vendored build
  reports a version whose hasher may not be the one it used is not answered, only accepted, and the
  consumer must present the result as an inference about the build.
- **It lists exact release tags, not version ranges.** A tag that is not listed returns no answer, so
  a dev or snapshot build fails to resolve rather than resolving to a neighbour. String equality is
  enough; the consumer needs no version-comparison rules.
- **It admits the gap.** Anything before 25.10.0 is unidentifiable and resolves to nothing at all.
  This replaces the `null` → THV_0 entry sketched earlier in this ADR, which would have claimed one
  answer for all of history and been wrong for most of it.

The table is generated, not written by hand. For each of the six boundary commits above, the releases
that contain it are computed with `git tag --contains`, so a change backported into an older
maintenance line produces a correct table rather than one inferred from version ordering. Run today,
no boundary commit has been backported, and every release maps to exactly one THV:

| THV | Releases |
| --- | --- |
| `std/v1.1` | 25.10.0 … 25.10.8 |
| `std/v1.2` | 25.11.0-edge, 25.12.0-edge |
| `std/v1.3` | 26.01.1-edge, 26.02.0-edge |
| `std/v1.4` | 26.03.0-edge … 26.03.4-edge, 26.04.0, 26.04.1 |
| `std/v1.5` | 26.04.2 … 26.04.6, 26.05.0-edge … 26.07.0-edge |
| `std/v1.6` | 26.08.0-edge |
| `std/v1.7` | 26.09.0-edge, 26.09.1-edge |

The generator must run against a full tag fetch. A release it does not see is simply absent, and an
absent release resolves to nothing, which is the correct failure.

### 7. What this does not cover

- **`HashMode`.** It is a runtime choice, not part of a THV. Two runs can share a THV and an
  identical field list and still hash differently because one ran under `DEEP` or `SHA256`. It
  remains a prerequisite below, unaddressed by this decision.
- **Patched and vendored builds.** They stamp correctly from now on; before the stamp they cannot be
  resolved honestly and must not be guessed.
- **Completeness of the diff.** Unchanged from the non-goals: the result explains, it does not prove.

### How the pieces land

In dependency order, each step standing on its own:

1. The `HashBuilder` encoding-rule fix. It is a bug on master independent of this ADR.
2. The spec files, `TaskHashSpec`, and the guard tests.
3. `TaskRun.hashVersion` and `AgentRun.hashVersion`, plus the two new lineage fields.
4. The `lineage` lists in the spec files.
5. The registry and its two sync guards.

## Prerequisites

The options above only describe the drift; they do not fix it. Three prerequisites were listed when
this ADR was a draft. The decision settles two of them and leaves one open.

**Settled by the decision:**

1. *Add the missing hash inputs to the lineage record.* Done with two fields, `resourcesBundle` and
   `envModules`. The stub-run marker needs no field of its own — a stub run already changes
   `codeChecksum`. The separate process-name field is dropped: `name` is a superset of the hashed
   process name, so the diff may over-report that key and can never under-report it.
2. *Ensure a new hash key cannot be added silently.* `StdSpecsGoldenTest` is that place. A key added
   to the hasher without a new spec makes the newest spec stop reproducing the default path, and the
   build fails. A key present in a spec but never populated fails the same test.

**Still open:**

3. **`HashMode` is not recorded.** It is a *runtime* choice, not part of a THV. Two runs can share a
   THV and an identical field list and still hash differently because one ran under `DEEP` or
   `SHA256`, and any mode added later widens that gap rather than closing it. What is needed is the
   *effective per-process* mode (`task.processor.getConfig().getHashMode()`), recorded as a field of
   the task run. The `mode` already present on every `Checksum` does not cover this: it is the global
   `NXF_CACHE_MODE` default, it describes the checksum it sits on rather than the task hash, and in
   records written before #7582 it is mislabelled as well. The cost is one enum-valued field per
   `TaskRun`, which is the cheapest item in this ADR and guards a whole class of otherwise invisible
   miss.

## Open questions

- **Lineage model version.** This work adds three fields to `TaskRun` (`hashVersion`,
  `resourcesBundle`, `envModules`) and one to `AgentRun`. That is read-compatible in both directions
  with Gson (new readers see `null`, old readers ignore unknown fields), and `LinTypeAdapterFactory`
  gates on strict equality of `version` (`lineage/v1beta1`), so no bump is strictly required. Whether
  growing the model without a bump is acceptable practice still needs settling, since `WorkflowRun`
  also feeds the execution hash (`CacheHelper.hasher(value).hash()` in `LinObserver.storeWorkflowRun()`)
  and new fields change that hash. Tracked by PR #7171.
- **Who generates `legacy-releases.json`, and when.** It must be regenerated on every release, from a
  full tag fetch, and a release missing from it resolves to nothing. Whether that belongs in the
  Nextflow release workflow or in a job owned by the schemas repository is not decided.
- **`-dump-hashes` alignment.** Whether the named-key output should be promoted to a documented,
  stable format consumable by tooling, or remain a debug aid. This design makes the two agree, since
  `-dump-hashes` prints the contributor names from the spec, which makes the question easier to answer
  either way.

Three questions from the draft are now answered:

- **THV form** — a semantic id, `std/<major>.<minor>`. The minor moves when the current code can
  still reproduce the older version, the major only when it cannot. The objection that "names can be
  rebound" is answered by the frozen fingerprint table: an edited spec changes its fingerprint and
  fails the build, so an id cannot quietly come to mean something else.
- **THV ownership across modules** — the encoding rules live in `nf-commons` as `EncodingRules`,
  next to the `HashBuilder` they configure; the key list lives in `nextflow`. A spec file names both,
  so one identifier covers both modules without either depending on the other.
- **Whether to wait on PR #6927** — no. The `TaskHasherFactory` extension point it needed landed on
  master independently (#7638), so this design builds on it directly.

## Links

- Refines [Data lineage](20250508-data-lineage.md)
- Related PR #6927 "Add versioned task hasher strategy" (`task-hasher-strategy`, open since March 2026,
  **not assumed to land**) — would introduce `TaskHasherFactory.Version` (`std/v1..v3`) and
  `NXF_TASK_HASH_VER`, and carries the module resources bundle fingerprint into the hash (#6914, closes
  #6128)
- Merged: #7638 (`631ed9734`) "Add extension points for task caching and object-storage access" —
  supplies the `TaskHasherFactory` extension point a versioned hasher needs, independently of #6927
- Proposed registry home: `nextflow-io/schemas`, under `task-hash/v1/`
- Related PR #7171 "ADR: lineage record version compatibility" (open) — covers the lineage model
  version question raised in the open questions above
- Fixed upstream: #7582 (`2c68fa1c9`) "Fix lineage checksum computed in standard mode regardless of
  label" — `Checksum.ofNextflow()` now hashes under the mode it labels
- Code: `modules/nextflow/src/main/groovy/nextflow/processor/TaskHasher.groovy`,
  `modules/nf-lineage/src/main/nextflow/lineage/model/v1beta1/TaskRun.groovy`,
  `modules/nf-lineage/src/main/nextflow/lineage/LinObserver.groovy`
