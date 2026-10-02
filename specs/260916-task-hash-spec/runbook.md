<!-- specs/260916-task-hash-spec/runbook.md -->
# Byte-exactness runbook

A cache hit on `-resume` is byte-exact hash equality, proven through `CacheDB` and
`TaskProcessor` rather than a unit assertion.

Run everything from `specs/260916-task-hash-spec/pipeline/`. `$PROTO` is the
prototype launcher (`../../../launch.sh`).

**Verification status of this document:** every section below was executed on
2026-10-02 against genuine released Nextflow distributions, and the quoted output is
the actual output. Docker and conda were both available, so the runs used the full
`-profile withtooling`: `P_CONTAINER` exercised the `CONTAINER` key against a real
image and `P_CONDA` the `CONDA` key against a resolved environment. Note conda lives
at `~/miniconda3/bin` and is not on `PATH` in a non-interactive shell — prepend it
before running, or `P_CONDA` silently runs without the conda key ever firing.

Keys exercised end to end: SESSION_ID, PROCESS_NAME, TASK_SOURCE, INPUTS,
EVAL_OUTPUTS, SCRIPT_VARS, BIN_ENTRIES, MODULE_BUNDLE, CONTAINER, CONDA, STUB_MARKER.
Still unexercised by the pipeline, by design: ENV_MODULES and SPACK (need HPC/spack
infrastructure); both are covered by unit oracles instead.

## 1. Positive: each spec reproduces its release

> **Check the baseline is a real release before trusting it.** Run
> `NXF_VER=<version> nextflow -version` and look at the build number. A genuine
> release has a real one (26.04.6 = `build 12646`); a locally installed snapshot
> reports `build 0`. Three distributions on this machine were such snapshots and had
> to be moved aside and re-downloaded: `26.08.0-edge`, `26.02.0-edge` and `25.10.0`.
> A baseline like that produces an apparent spec failure that is really a mislabelled
> binary. Delete or ignore such distributions when validating.



| Spec | Baseline | Result (2026-10-02) |
| --- | --- | --- |
| `std/v1.1` | `NXF_VER=25.10.0` (build 10289) | 8/8 CACHED |
| `std/v1.2` | `NXF_VER=25.11.0-edge` (build 10558) | 8/8 CACHED |
| `std/v1.3` | `NXF_VER=26.02.0-edge` (build 11371) | 8/8 CACHED |
| `std/v1.4` | `NXF_VER=26.04.1` (build 12112) | 8/8 CACHED |
| `std/v1.5` | `NXF_VER=26.04.6` (build 12646) | 8/8 CACHED |
| `std/v1.6` | `NXF_VER=26.08.0-edge` (build 13213) | 8/8 CACHED |
| `std/v1.7` | `NXF_VER=26.09.0-edge` (build 13644) | 8/8 CACHED |
| default path (no env var) | `NXF_VER=26.09.0-edge` | 8/8 CACHED |

The last row is the default-unchanged guarantee: a run laid down by a released
Nextflow still resumes under the prototype with no spec selected.

`std/v1.1` additionally needs the asset-root fixture of section 1c — the main
pipeline cannot trigger the behaviour it restores.

For each row:

```bash
rm -rf work .nextflow*
NXF_VER=<version> nextflow run main.nf -profile withtooling
NXF_TASK_HASH_VER=<spec> $PROTO run main.nf -profile withtooling -resume
```

**Pass:** every process reports `CACHED` in the second run.

### Self-resume — verified

```
$ $PROTO run main.nf
[SUCCESS] completed=8 failed=0 cached=0
$ $PROTO run main.nf -resume
[SUCCESS] completed=0 failed=0 cached=8
```

## 1b. Cross-spec discrimination WITHOUT historical releases

Section 3 below needs old Nextflow distributions. This variant needs none, and was
run on 2026-09-22 against the prototype: lay down a `std/v1.7` baseline, then resume
under each older spec. Each step back must miss exactly the key the corresponding
upstream commit changed.

```bash
for spec in std/v1.6 std/v1.5 std/v1.3; do
  rm -rf work .nextflow*
  $PROTO run main.nf -profile withtooling                       # std/v1.7 baseline
  NXF_TASK_HASH_VER=$spec $PROTO run main.nf -profile withtooling -resume
done
```

**Pass:** a monotone ladder — observed exactly:

| resumed under | processes that missed | maps to |
| --- | --- | --- |
| `std/v1.6` | `P_EVAL` | #7575, eval-output form |
| `std/v1.5` | + `P_MODULE_BUNDLE` | #6914, module-bundle key |
| `std/v1.3` | + `P_MAP_INPUT` | #6679, `orderIndependentMaps` |

This proves those specs are mutually distinct in exactly the predicted places. It
does NOT replace section 1: only a run against a real released Nextflow shows a spec
reproduces that *release* rather than merely differing correctly from its siblings.

## 1c. Asset-root fixture — the only way to test `std/v1.1`

`std/v1.1` restores the pre-#6605 `isAssetFile`, which looked at `baseDir` alone.
The main pipeline cannot trigger it: every file it hashes is already under `baseDir`.
The difference only appears for a file inside the project repository but outside
`baseDir`, which needs the main script to live in a subdirectory.

Build the fixture once:

```bash
R=~/.nextflow/assets/testorg/hashv11
mkdir -p $R/sub $R/shared
echo "manifest.mainScript = 'sub/main.nf'" > $R/nextflow.config
echo "asset payload" > $R/shared/data.txt
cat > $R/sub/main.nf <<'NF'
process P_ASSET {
  input:
  path x
  output:
  stdout
  script:
  "cat $x"
}

workflow {
  P_ASSET( file("${projectDir}/../shared/data.txt") ) | view
}
NF
( cd $R && git init -q && git add -A && git -c user.email=t@t -c user.name=t commit -qm init \
  && git remote add origin https://github.com/testorg/hashv11.git )
```

A git repository and a remote are both required: `isAssetFile` returns false when
`session.commitId` is null. Run everything with `NXF_OFFLINE=true` so Nextflow does
not try to reach the fake remote.

```bash
rm -rf work .nextflow*
NXF_OFFLINE=true NXF_VER=<baseline> nextflow run testorg/hashv11
NXF_OFFLINE=true NXF_TASK_HASH_VER=<spec> $PROTO run testorg/hashv11 -resume
```

**Pass:** a two-way split — each spec reproduces its own era and misses the other.
Observed 2026-10-02:

| Baseline | `std/v1.1` | `std/v1.2` |
| --- | --- | --- |
| 25.10.0 | CACHED | re-runs |
| 25.11.0-edge | re-runs | CACHED |
| 25.12.0-edge | re-runs | CACHED |
| 26.01.1-edge | re-runs | CACHED |

A run where both specs cache means the flag is not reaching the file — see the
`nested()` propagation in `HashBuilder`. Confirm by touching `shared/data.txt`:
`std/v1.1` hashes metadata so the digest must move, `std/v1.2` hashes content so it
must not.

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
process stays `CACHED`. All three were run against std/v1.7 self-resume and matched
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
NXF_VER=26.08.0-edge nextflow run main.nf -profile withtooling
NXF_TASK_HASH_VER=std/v1.7 $PROTO run main.nf -profile withtooling -resume
```

**Pass:** everything `CACHED` **except `P_EVAL`**, which re-runs.

```bash
rm -rf work .nextflow*
NXF_VER=26.04.6 nextflow run main.nf -profile withtooling
NXF_TASK_HASH_VER=std/v1.6 $PROTO run main.nf -profile withtooling -resume
```

**Pass:** only `P_MODULE_BUNDLE` re-runs.

```bash
rm -rf work .nextflow*
NXF_VER=26.02.0-edge nextflow run main.nf -profile withtooling
NXF_TASK_HASH_VER=std/v1.5 $PROTO run main.nf -profile withtooling -resume
```

**Pass:** `P_MAP_INPUT` re-runs. If *nothing* re-runs, the Map never reached the
hasher and the encoding dimension is untested — fix the fixture before trusting
`std/v1.3`.

## 4. Stub marker

```bash
rm -rf work .nextflow*
NXF_VER=26.08.0-edge nextflow run main.nf -profile withtooling -stub-run
NXF_TASK_HASH_VER=std/v1.6 $PROTO run main.nf -profile withtooling -stub-run -resume
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
[P_BASIC] cache hash: 100ecab06f46de7badf101d2818a692e; spec: std/v1.7; entries: [
    { "spec": "std/v1.7", "fingerprint": "7bbb059fc53596d9e06281d53275d3ed" },
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
