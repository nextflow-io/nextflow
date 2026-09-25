<!-- specs/260916-task-hash-spec/runbook.md -->
# Byte-exactness runbook

A cache hit on `-resume` is byte-exact hash equality, proven through `CacheDB` and
`TaskProcessor` rather than a unit assertion.

Run everything from `specs/260916-task-hash-spec/pipeline/`. `$PROTO` is the
prototype launcher (`../../../launch.sh`).

**Verification status of this document:** sections 1 and 3 depend on historical
Nextflow releases (`NXF_VER=26.02.0-edge` etc.) that download old distributions —
an agent cannot fetch those, so those rows are written for a human to run and are
marked accordingly below. Everything that does not require a historical release —
the `std/v4` self-resume in section 1, all three negative controls in section 2,
and the stub marker in section 4 — was run against this build with `$PROTO` and
its actual output is quoted inline. Docker and conda were BOTH available on re-verification
(2026-09-22), so the runs below used the full `-profile withtooling`: `P_CONTAINER`
exercised the `CONTAINER` key against a real image and `P_CONDA` the `CONDA` key
against a resolved environment. Note conda lives at `~/miniconda3/bin` and is not on
`PATH` in a non-interactive shell — prepend it before running, or `P_CONDA` silently
runs without the conda key ever firing.

Keys exercised end to end: SESSION_ID, PROCESS_NAME, TASK_SOURCE, INPUTS,
EVAL_OUTPUTS, SCRIPT_VARS, BIN_ENTRIES, MODULE_BUNDLE, CONTAINER, CONDA, STUB_MARKER.
Still unexercised by the pipeline, by design: ENV_MODULES and SPACK (need HPC/spack
infrastructure); both are covered by unit oracles instead.

## 1. Positive: each spec reproduces its release

> **Check the baseline is a real release before trusting it.** Run
> `NXF_VER=<version> nextflow -version` and look at the build number. A genuine
> release has a real one (26.04.6 = `build 12646`); a locally installed snapshot
> reports `build 0`. On this machine `~/.nextflow/framework/26.08.0-edge` and
> `26.02.0-edge` are both `build 0` and contain code that does NOT match the tags
> of those names — the "26.08.0-edge" binary hashes eval outputs as a raw map,
> which is post-#7575 behaviour, while the v26.08.0-edge tag predates #7575. A
> baseline like that produces an apparent spec failure that is really a mislabelled
> binary. Delete or ignore such distributions when validating.



| Spec | Baseline | Verified here? |
| --- | --- | --- |
| `std/v1` | `NXF_VER=26.02.0-edge` — 7/8, `P_BASIC` misses by design, see note below | No — human run required |
| `std/v2` | `NXF_VER=26.04.6` — VERIFIED 2026-09-22, 8/8 CACHED | No — human run required |
| `std/v3` | `NXF_VER=26.08.0-edge` — VERIFIED 2026-09-22, 8/8 CACHED (build 13213) | No — human run required |
| `std/v4` | the prototype itself, no env var | **Yes** — see below |

For each row:

```bash
rm -rf work .nextflow*
NXF_VER=<version> nextflow run main.nf -profile withtooling
NXF_TASK_HASH_VER=<spec> $PROTO run main.nf -profile withtooling -resume
```

**Pass:** every process reports `CACHED` in the second run.

### std/v4 row — verified

```
$ $PROTO run main.nf
[SUCCESS] completed=8 failed=0 cached=0
$ $PROTO run main.nf -resume
[SUCCESS] completed=0 failed=0 cached=8
```

All 8 processes cached on the second run with no env var (std/v4 is the default
spec) — this is the default-unchanged guarantee.


> **`std/v1` expected result is 7/8, not 8/8.** `P_BASIC` re-executes against
> `26.02.0-edge`, and this is correct behaviour rather than a spec defect. Commit
> `785e801ad` (#7165, 2026-05-21) made the strict parser fold referenced `params.*`
> into the task global vars, so `SCRIPT_VARS` today hashes
> `[params.greeting=hello, task.ext.flavour=vanilla]` where February hashed
> `[task.ext.flavour=vanilla]`. Every other entry of `P_BASIC` is byte-identical.
> A spec fixes which keys are hashed and how they are encoded; it cannot restore how
> an era *derived* a value, because contributors call today's helpers. Treat a
> `P_BASIC`-only miss under `std/v1` as a pass; any OTHER process missing is a real
> failure.

## 1b. Cross-spec discrimination WITHOUT historical releases

Section 3 below needs old Nextflow distributions. This variant needs none, and was
run on 2026-09-22 against the prototype: lay down a `std/v4` baseline, then resume
under each older spec. Each step back must miss exactly the key the corresponding
upstream commit changed.

```bash
for spec in std/v3 std/v2 std/v1; do
  rm -rf work .nextflow*
  $PROTO run main.nf -profile withtooling                       # std/v4 baseline
  NXF_TASK_HASH_VER=$spec $PROTO run main.nf -profile withtooling -resume
done
```

**Pass:** a monotone ladder — observed exactly:

| resumed under | processes that missed | maps to |
| --- | --- | --- |
| `std/v3` | `P_EVAL` | #7575, eval-output form |
| `std/v2` | + `P_MODULE_BUNDLE` | #6914, module-bundle key |
| `std/v1` | + `P_MAP_INPUT` | #6679, `orderIndependentMaps` |

This proves the four specs are mutually distinct in exactly the predicted places. It
does NOT replace section 1: only a run against a real released Nextflow shows a spec
reproduces that *release* rather than merely differing correctly from its siblings.

## 2. Negative controls

Without these, a run where everything caches for unrelated reasons reads as a pass.
After a passing step 1, each of these must make exactly the named process re-execute:

```bash
# change a script body
sed -i 's/echo eval-process/echo eval-process-CHANGED/' main.nf   # P_EVAL re-runs
# change a value input
sed -i "s/Channel.value('word')/Channel.value('WORD')/" main.nf   # P_BASIC re-runs
# rename a process
sed -i 's/P_BIN/P_BIN_RENAMED/g' main.nf                          # P_BIN re-runs
```

Revert each change before the next.

**Pass:** the named process reports as executed (not `CACHED`) and every other
process stays `CACHED`. All three were run against std/v4 self-resume and matched
exactly:

```
# eval body change
[PROCESS ...] P_EVAL
[SUCCESS] completed=1 failed=0 cached=7

# value input change
[PROCESS ...] P_BASIC
[SUCCESS] completed=1 failed=0 cached=7

# process rename
[PROCESS ...] P_BIN_RENAMED
[SUCCESS] completed=1 failed=0 cached=7
```

Each change reverted afterwards; a clean `-resume` returns to `cached=8`.

## 3. Cross-spec discrimination

The strongest check: the specs must differ in exactly the predicted place. All
three rows below require a historical `NXF_VER` and are written for a human to
run; they were not executed by the agent.

**On `P_BAG_INPUT` — verified and dropped.** The plan called for a second H1-vs-H2
discriminator: a multi-file `path` input, delivered as an `ArrayBag<FileHolder>`,
to independently exercise `cacheFunnelFirst` (the other half of `EncodingRules`,
alongside `orderIndependentMaps`) on the premise that `ArrayBag` implements both
`Bag` and `CacheFunnel`. Checked against
`modules/nextflow/src/main/groovy/nextflow/util/ArrayBag.groovy`: it implements
`Bag<E>, List<E>, KryoSerializable` — **not** `CacheFunnel`. A further check across
the codebase (`FileHolder`, `GroupKey`, `SecretImpl`, `CmdLineOptionMap`,
`PluginRef`, `AzPoolOpts`/`AzFileShareOpts` — every `CacheFunnel` implementor) found
no type that is simultaneously a `Map`/`Bag`/`Set` and a `CacheFunnel`. In
`HashBuilder.with()` (`modules/nf-commons/src/main/nextflow/util/HashBuilder.java`),
`cacheFunnelFirst` only changes the outcome for such a dual object — against a bag
of plain `FileHolder`s the funnel-vs-collection branch order is unobservable.
Separately, `StdSpecs` never varies `cacheFunnelFirst` independently of
`orderIndependentMaps`: `EncodingRules.LEGACY` sets both `false` and
`RECORD_TYPES` sets both `true`, so no pair of the four shipped specs could isolate
`cacheFunnelFirst` even given a dual object. A `P_BAG_INPUT` process was not added
to `main.nf` — it would not have discriminated anything, and a step that cannot
fail is not a test. Net effect: `P_MAP_INPUT` is the only available H1-vs-H2
discriminator, and it covers `orderIndependentMaps` only. `cacheFunnelFirst` is
untested by this pipeline; see "Known limitation" below.

```bash
rm -rf work .nextflow*
NXF_VER=26.08.0-edge nextflow run main.nf -profile withtooling      # H3
NXF_TASK_HASH_VER=std/v4 $PROTO run main.nf -profile withtooling -resume
```

**Pass:** everything `CACHED` **except `P_EVAL`**, which re-runs.

```bash
rm -rf work .nextflow*
NXF_VER=26.07.0-edge nextflow run main.nf -profile withtooling      # H2
NXF_TASK_HASH_VER=std/v3 $PROTO run main.nf -profile withtooling -resume
```

**Pass:** only `P_MODULE_BUNDLE` re-runs.

```bash
rm -rf work .nextflow*
NXF_VER=26.02.0-edge nextflow run main.nf -profile withtooling      # H1
NXF_TASK_HASH_VER=std/v2 $PROTO run main.nf -profile withtooling -resume
```

**Pass:** `P_MAP_INPUT` re-runs. If *nothing* re-runs, the Map never reached the
hasher and the encoding dimension is untested — fix the fixture before trusting
`std/v1`.

## 4. Stub marker

```bash
rm -rf work .nextflow*
NXF_VER=26.08.0-edge nextflow run main.nf -profile withtooling -stub-run
NXF_TASK_HASH_VER=std/v3 $PROTO run main.nf -profile withtooling -stub-run -resume
```

**Pass:** `P_STUB` reports `CACHED`.

### Same-build stub self-resume — verified

The historical-baseline row above needs a human; the same-build case does not,
and was run:

```
$ $PROTO run main.nf -stub-run
[SUCCESS] completed=8 failed=0 cached=0
$ $PROTO run main.nf -stub-run -resume
[SUCCESS] completed=0 failed=0 cached=8
```

`P_STUB` (and everything else) reports `CACHED`, confirming `hasStubBlock()` /
`STUB_MARKER` round-trips through a real resume.

## When resume fails

Resume can fail for reasons unrelated to hashing — a February build's LevelDB/Kryo
cache may not be readable by a master build. Use the independent second track to
tell the two apart:

```bash
NXF_VER=<version> nextflow run main.nf -profile withtooling -dump-hashes json 2>&1 | tee baseline.log
NXF_TASK_HASH_VER=<spec> $PROTO run main.nf -profile withtooling -dump-hashes json 2>&1 | tee proto.log
```

Compare the reported `cache hash:` per process. Equal hashes with a failing resume
means the cache plumbing, not the spec. Unequal hashes: the prototype's per-key
entries name the key that differs.

`-dump-hashes json` was confirmed to run and to emit exactly this shape (spec id,
fingerprint, and one entry per contributing key) against the current build:

```
$ $PROTO run main.nf -dump-hashes json
[P_BASIC] cache hash: 100ecab06f46de7badf101d2818a692e; spec: std/v4; entries: [
    { "spec": "std/v4", "fingerprint": "7bbb059fc53596d9e06281d53275d3ed" },
    { "key": "SESSION_ID", "hash": "a0907d420d6a705c6278c11daed8d36a" },
    { "key": "PROCESS_NAME", "hash": "3cd9d61e7555804d3902059aa1d57e78" },
    { "key": "TASK_SOURCE", "hash": "6341d19c030e8f7796cf7d8416f09b4c" },
    { "key": "INPUTS", "hash": "9880fc64a65fa763df4c2e732fa16339" },
    { "key": "SCRIPT_VARS", "hash": "58f6bb9b99bd0e5da33e85a02a3c37b3" }
]
```

## Container digest

`P_CONTAINER` is pinned to
`quay.io/nextflow/bash@sha256:bea0e244b7c5367b2b0de687e7d28f692013aa18970941c7dd184450125163ac`.
This digest was resolved and confirmed to exist by pulling `quay.io/nextflow/bash:latest`
and reading back its `RepoDigests` entry — it is real, not a placeholder. If it
stops resolving by the time this runbook is run, re-pull the image and substitute
the current digest; keep it pinned by digest, or the `CONTAINER` key varies for
environmental reasons unrelated to the spec under test.

## Known limitation

The era map is derived from the history of `TaskHasher` and `HashBuilder` only. A
change in a `CacheFunnel` implementor (`FileHolder`, `ArrayBag`, …) would move hashes
without touching either file. If a positive run misses, either the spec is wrong or
there is a boundary we have not found — `-dump-hashes` says which key.

Additionally (see section 3): `cacheFunnelFirst` has no pipeline fixture that
discriminates it independently of `orderIndependentMaps`, because no type in the
current codebase is both a `Map`/`Bag`/`Set` and a `CacheFunnel`, and the four
shipped specs never vary the two flags independently. This is a real gap in what
this runbook can prove, not an oversight in the fixture — closing it would require
either a new dual-purpose type in production code or a spec that flips the flags
independently, both out of scope for this validation pass.
