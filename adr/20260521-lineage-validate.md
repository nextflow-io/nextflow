# Lineage Validate: Semantic Equivalence Check for Workflow Runs

- Authors: Edmund Miller
- Status: Draft
- Deciders: TBD
- Date: 2026-05-21
- Tags: lineage, cli, reproducibility, ci, testing

Technical Story: Surface a `nextflow lineage validate` command that lets users assert two workflow runs are semantically equivalent — for CI regression checks, reproducibility audits, and divergence debugging.

## Summary

`nextflow lineage validate` compares two workflow runs (or a run against a saved snapshot) and reports whether they are *semantically equivalent* — same key inputs, same published outputs (by checksum), same workflow identity. The command is the CLI face of a shared `LineageValidator` core also used by `LineageSnapshotter` (Spock integration). Comparison is **outputs-flat with one-hop causality enrichment**; runtime fingerprint (containers, executor) is informational by default. Platform (Seqera) baselines are resolved via a pluggable SPI implemented by `nf-tower`.

## Problem Statement

Pipelines that depend on data lineage need a way to assert "two runs produced the same result." Today the lineage store records workflow runs, task runs, file outputs, parameters, and workflow identity (repository + commit), but there is no canonical way to ask "did anything that matters change?" Users currently must:

- Inspect lineage JSON manually,
- Write ad-hoc Groovy/Python scripts (e.g., Stimulus' `seqera-compare-runs` skill),
- Or trust pipeline-level smoke tests that don't see the data-layer differences.

This blocks three converging use cases:

1. **CI regression checks** — after a refactor or version bump, assert outputs and key inputs didn't drift.
2. **Reproducibility audits** — certify that a published run matches a re-run.
3. **Divergence debugging** — "why is run B different from run A?"

The first ships of `nextflow lineage validate` (commits 70b62f2e0, 7e2239a6b) and `LinNormalizer` / `LineageSnapshotter` cover the basic two-run-by-LID path. This ADR pins the design before that surface stabilises — defining the equivalence model, snapshot file format, CI ergonomics, and Platform integration seam.

## Goals or Decision Drivers

- **CI-first defaults** — the dominant use case. Fast pass/fail, scriptable exit codes, structured `--json` for tooling, auto-detection of CI environments.
- **Audit and debugging via flags** — same command, opt-in flags surface strictness (`--strict-runtime`, `--depth`) and drill-down (`--trace`).
- **Single source of truth** — the CLI, the Spock snapshotter, and any future Platform integration share one `LineageValidator` core. No drift between "what the test asserts" and "what CI checks."
- **Pluggable baselines** — accept `lid://`, a local snapshot file, or a Platform run (`tower://`) without coupling `nf-lineage` to Platform auth.
- **Honest semantics** — equivalence = recorded checksums + key inputs + workflow identity. Validate does not re-hash files, prove repeatability, or build a full causality graph.
- **Stable surface for pipeline authors** — `--ignore-fields` paths and the snapshot format are versioned and documented.

## Non-goals

- **Bit-for-bit content checking.** Validate trusts recorded checksums; it does not re-hash or open files. If a pipeline records a wrong checksum, validate cannot catch it.
- **A reproducibility guarantee.** Validate asserts that two *recorded* lineage trees are equivalent. It does not prove the runs would re-execute identically (reproducibility = repeatability + lineage match; validate covers only the second half).
- **A full DAG / causality engine.** Comparison flattens to published outputs + key inputs; the report enriches diverging outputs with one-hop upstream context. Deep DAG analysis is a separate `lineage explain` command (future work).
- **Platform-aware by itself.** Core `nextflow lineage validate` is local-only. `tower://` and other remote references resolve only when a resolver plugin is installed. Without `nf-tower`, `tower://` refs fail with a clear error.
- **Multi-baseline (N-to-1) canary validation in v1.** Single baseline only; quorum semantics are deferred.
- **Runtime reproducibility by default.** Container digests, conda hashes, executor types are *informational* in the diff. Strict assertion requires `--strict-runtime`.
- **Long-term archival snapshot format.** Snapshot files are schema-versioned and migratable across minor Nextflow versions but are not promised to be readable across major releases without migration.

## Design Decisions

The grilling that produced this ADR walked 16 design questions; below are the decisions and their rationale.

### D1 — Primary use case: CI-led, three-use-case design

Three converging use cases (CI, audit, debug) are first-class. **CI is the lead**: defaults optimise for fast pass/fail, structured output, sensible behaviour without flags. Audit and debug are served by `--strict-*`, `--trace`, and `--depth` flags.

### D2 — Equivalence unit: outputs + key inputs

"Two runs are equivalent" means: same key inputs (params, repository, commit) produced the same published outputs (by recorded checksum). The internal task graph is implementation detail. Refactors that don't change outputs don't fail validation. This aligns with what the lineage model already represents as the user-facing contract.

### D3 — File equivalence: checksum-only

Each `FileOutput` carries a checksum (algorithm, provider, value) in the lineage record. Validate compares checksums. It does **not**:

- Open files to compute fresh hashes
- Run type-aware (VCF/BAM/text) diffs
- Re-hash with a different algorithm

Non-determinism in outputs (embedded timestamps, parallel-write ordering, compression metadata) is a pipeline-design problem; the pipeline author is responsible for producing stable checksums (e.g., `gzip -n`, deterministic sort). Validate exposes `--ignore-fields outputs.<name>.checksum` as the escape hatch for legitimately volatile outputs.

### D4 — Baseline model: peer LID *and* snapshot file from day one

`--against` accepts either an LID (`lid://wf-abc123`) or a path to a snapshot file. A complementary `--save-snapshot <path>` writes the normalized tree of the *new* run for use as a future baseline. The snapshot format aligns with what `LineageSnapshotter` already writes.

Future work: tagged/named baselines (`lid://golden/v1.0` aliases) and multi-baseline quorum.

### D5 — Subworkflow handling: flatten to terminal outputs

Validate compares the *top-level* workflow run's published outputs. Sub-workflow structure is an implementation detail — if the published outputs match, the run is equivalent. This avoids spurious failures when refactors move logic between sub-workflows but keep outputs stable.

Aligns with D2 (outputs + key inputs).

### D6 — Failure output: human diff default, `--json` for CI

- **Default**: JGit-style unified diff to stdout, grouped by category (see D13). Matches current behaviour.
- **`--json`**: stable structured schema. Array of differences; each carries `category`, `path`, `expected`, `actual`, optionally `cause` (see D15).
- **`--summary`**: counts-only one-liner (`5 outputs differ, 1 param differs`).
- **Exit codes**: `0` = equivalent, `1` = differs, `2` = error (load failure, schema mismatch, resolver error).

### D7 — Auto-detect CI environments

Validate inspects environment variables to choose a default output mode:

| Env | Default mode |
|-----|--------------|
| `CI=true` or `GITHUB_ACTIONS=true` | `--json` to stdout |
| `AGENT=1` or `CLAUDECODE=1` | one-line summary + `--json` tail |
| `GITHUB_ACTIONS=true` (additionally) | emit `::error::` / `::notice::` annotations and write JSON to `$GITHUB_STEP_SUMMARY` if set |
| (none) | human-readable diff |

Explicit `--format human|json|summary` always overrides. Precedence is documented and stable.

### D8 — Ignore mechanism: flag + config, JSONPath-style

- **Flag**: `--ignore-fields params.outdir,outputs.*.modifiedAt,workflow.commitId`
- **Config**: `nextflow.config` block:
  ```groovy
  lineage {
      validate {
          ignore = ['params.outdir', 'outputs.*.modifiedAt']
      }
  }
  ```
- **Built-in `EPHEMERAL_FIELDS`** always apply (sessionId, createdAt, modifiedAt, name, source, workflowRun, taskRun, path), unaffected by user config.
- **Syntax**: dot-paths with `*` (single segment) and `**` (recursive) wildcards. No JSONPath expression engine — keep matching simple and predictable.

### D9 — Workflow identity: `commitId` + `repository`

Two runs are "the same pipeline" if they share repository URL and commit ID. Script file checksums are not compared when commit ID is present (Git invariant guarantees content match). **Dirty-tree fallback**: when no commit ID is recorded (running from a working tree), validate falls back to comparing `scriptFiles` checksums.

Container image digests / conda env hashes are *runtime identity*, addressed in D17.

### D10 — Schema evolution: refuse on major mismatch, normalize on minor

Snapshots and lineage records carry `schemaVersion` (currently `v1beta1`). Comparison rules:

- **Identical** schema → direct compare.
- **Same major, different minor** → normalize both sides; fields present only on the newer side are dropped (treated as ignored).
- **Different major** (e.g., `v1` vs `v2`) → hard error, exit code 2, message naming the version mismatch and the recommended migration path.

This sets up the lineage model's `v1beta*` lifecycle: minor bumps add fields, major bumps may rename or remove and require migration.

### D11 — CLI vs `LineageSnapshotter`: extract shared `LineageValidator`

Both the CLI command (`LinCommandImpl.validate`) and `LineageSnapshotter` (Spock integration) call a single `LineageValidator` core:

```
LineageValidator
  ├─ resolve(ref) → LineageTree     // via LineageResolver SPI
  ├─ normalize(tree, options)
  ├─ compare(treeA, treeB) → DifferenceReport
  └─ render(report, format)         // human | json | summary | gh-actions
```

The CLI is a thin wrapper that picks up flags / env / config. `LineageSnapshotter` is a thin wrapper that knows about Spock conventions and `UPDATE_SNAPSHOTS`. The comparison algorithm and snapshot format live in `LineageValidator` only.

This avoids the failure mode where "the test passed but CI failed" (or vice versa).

### D12 — Platform integration: pluggable `LineageResolver` SPI

`LineageValidator.resolve(ref)` dispatches on the URI scheme:

- `lid://...` → built-in `LinStore` resolver
- file path / `file://...` → built-in `SnapshotFileResolver`
- `tower://<runId>` → contributed by `nf-tower` plugin (joins Platform API `/workflow/:id` + `/workflow/:id/tasks` against the local lineage store on task hash)
- `stimulus://<fixture>` → potentially contributed by a future test fixture plugin

Without `nf-tower` installed, `tower://` references fail with: `tower:// references require the nf-tower plugin (install with: nextflow plugin install nf-tower)`.

The resolver SPI is the design seam Stimulus' `seqera-compare-runs` skill validates: task-hash join (primary) + workdir fallback for matching tasks across runs that have different LIDs.

### D13 — Difference categories: outputs / params / workflow-identity / resources

Each difference is bucketed:

| Category | Examples | Fails validation? |
|---|---|---|
| `outputs` | FileOutput checksum, missing/extra output | Yes |
| `params` | Workflow parameter value | Yes |
| `workflow-identity` | Repository, commit, scriptFiles (dirty tree) | Yes |
| `resources` | Peak memory, CPU time, wall time | No (informational) |
| `runtime` | Container digest, conda hash, executor | No (informational; see D17) |
| `status` | Workflow run status (SUCCESS / FAILED) | Yes (see D16) |

Categories appear as JSON arrays and as section headers in human output. The split mirrors Stimulus' `status / hash / resource / lineage` split.

### D14 — Output join key: relative path under `-outputDir`

Two `FileOutput` records are "the same output" if they share the same path *relative to* the workflow run's `-outputDir` (recorded in the lineage). This is stable across runs in different working directories and matches the user's mental model: "the published outputs at the same relative paths should be equivalent."

If `-outputDir` is not in the lineage record (older runs), fall back to the longest common path suffix — best-effort with a warning.

### D15 — Causality in failure reports: one-hop upstream

For each diverging `outputs` entry, the report includes:

- The producing `TaskRun` LID and content hash
- The diverging input(s) that caused the task to produce a different output (if findable in one upstream hop)

JSON schema:
```json
{
  "category": "outputs",
  "path": "results/sample1/aligned.bam",
  "expected": { "checksum": "sha256:abc..." },
  "actual":   { "checksum": "sha256:def..." },
  "cause": {
    "taskRun": "lid://t-xyz",
    "taskHash": "12/abcdef",
    "differingInputs": [ "params.reference" ]
  }
}
```

Comparison itself stays output-flat (CI-fast); the *report* does the one-hop walk for diverging entries only. Full upstream traversal is `--trace` (debug mode) or a future `lineage explain` command.

### D16 — Failed runs: compare what exists, surface status divergence

Validate does not refuse to compare a failed run against a successful one. Stimulus' failure-isolation use case explicitly needs this. Behaviour:

- Workflow run status (SUCCESS / FAILED / etc.) is a comparison field in the `status` category — divergence is reported and fails validation.
- For each tree, walk whatever lineage exists. Missing task runs / outputs on one side are reported in their respective categories (`outputs: missing on B`).
- One-hop causality (D15) is best-effort when task records exist.

### D17 — Runtime identity: informational, opt-in strictness

Container image digests, conda environment hashes, executor types — when recorded on `TaskRun` — appear in the `runtime` category. By default, runtime drift is reported but does not fail validation (CI ergonomics: container patches happen frequently). `--strict-runtime` promotes the `runtime` category to a fail-if-different category.

If the lineage model does not yet capture a particular runtime field, validate reports it as `unknown` rather than asserting equivalence. The lineage model expansion is a separate work item (see "Model Gaps" below).

### D18 — Non-file outputs: model gap, scope honestly

The current `v1beta1` lineage model captures `FileOutput` but not value emits (primitive return values, value channels, JSON-like emits). For v1 of validate:

- Validate compares whatever the model records — `FileOutput`s, `Parameter`s in `WorkflowRun`, etc.
- Value emits are **out of scope** until the lineage model adds `ValueOutput` (named follow-up).
- The ADR explicitly enumerates what's covered and what's not; pipeline authors who care about a value emit can persist it as a file or a parameter.

## Considered Options (for the core design seam)

### Option A — Direct port of the current implementation

Keep `LinCommandImpl.validate` and `LineageSnapshotter` as parallel implementations sharing only `LinNormalizer`. Add `--json` and snapshot-file support to each.

- Good, because lowest churn — extends what exists.
- Bad, because the comparison and snapshot logic drift between the two over time.
- Bad, because adding a Platform resolver requires plumbing twice.
- Bad, because the test-vs-CI parity failure mode is real and silent.

**Rejected** — drift cost outweighs short-term simplicity.

### Option B — `LineageValidator` core + resolver SPI (chosen)

Extract a shared `LineageValidator` with a `LineageResolver` SPI. CLI and `LineageSnapshotter` become thin wrappers. `nf-tower` contributes a `tower://` resolver.

- Good, because one place defines equivalence and snapshot format.
- Good, because Platform integration is a plugin concern, not a core concern.
- Good, because matches the existing Nextflow plugin architecture (`nf-amazon`, `nf-google`, `nf-tower` all add CLI/resolver hooks).
- Good, because Stimulus' `seqera-compare-runs` skill can be re-modelled around the SPI.
- Bad, because more surface area in `nf-lineage` (validator interface, resolver SPI, snapshot format spec).
- Bad, because plugin developers must learn another extension point.

**Chosen.**

### Option C — Move validate entirely into `nf-tower` / a new plugin

Pull validate out of core. CLI stays as a stub that delegates to whichever validate plugin is installed.

- Good, because keeps core lean.
- Bad, because validate is a *first-class* lineage feature; gating it on a plugin install is bad UX.
- Bad, because `LineageSnapshotter` is in `nf-lineage` already; moving validate out splits the abstraction.
- Bad, because makes the CLI command's behaviour dependent on plugin availability.

**Rejected** — validate belongs in core; the *Platform* part belongs in a plugin.

## Solution / Decision Outcome

Ship `nextflow lineage validate` backed by a shared `LineageValidator` core in `nf-lineage`, with comparison flattened to published outputs + key inputs + one-hop causality. Defaults are CI-friendly (auto-detect, exit codes, structured `--json`). Platform run baselines (`tower://`) are resolved through a `LineageResolver` SPI implemented by `nf-tower`. Runtime fingerprints are informational unless `--strict-runtime`.

## Rationale & Discussion

### Why outputs-flat (and not full DAG)

A full-DAG comparison is technically more accurate but practically hostile: every harmless refactor (renamed process, reordered channel ops, sub-workflow restructuring) becomes a validate failure. Pipeline authors will start `--ignore-fields`-ing huge swaths of the lineage, which defeats the safety net. Outputs are the *contract* between a pipeline and its consumers; comparing the contract is what users actually want from a CI check.

For the audit use case that *does* need full-DAG comparison, a future `--depth all` or sibling `lineage diff --deep` can target that case. The two designs are compatible.

### Why one-hop causality

Output-flat comparison answers "did anything diverge?" One-hop causality answers "*why*?" for the divergent outputs only — which is what the human reading the failure report needs. Deeper walks belong in a dedicated explain/trace command, not in the default validate report (cost and noise grow nonlinearly with depth).

### Why a resolver SPI rather than baking Platform in

The Stimulus skill (seqeralabs/stimulus#535) already demonstrates the value of joining Platform runs with lineage on task hash. Hard-coding a Platform client in `nf-lineage` couples the core to a specific deployment of Seqera Platform, breaks for self-hosted Tower instances with different auth, and makes the OSS path harder to test in isolation. The SPI keeps `nf-lineage` self-contained and lets `nf-tower` (already a plugin) own Platform concerns.

### Why categories with informational tiers

A flat list of diffs forces a binary outcome on every difference. Container image digests change every week in many deployments; treating that as a hard fail makes validate noisy enough that teams disable it. Categorising (`outputs` is contract, `resources` and `runtime` are informational) lets CI fail on contract drift without being dragged down by infrastructure churn. Users who *want* strict runtime checks opt in.

### Why snapshot format = normalized tree + version

Storing the normalized tree (rather than raw records) means snapshots stay readable even when the lineage record format adds fields. Storing a `schemaVersion` makes it possible to migrate (or refuse) when the model evolves. Storing the original LID, Nextflow version, and timestamp gives auditors traceability. This shape is exactly what `LineageSnapshotter` already writes; the ADR formalises it.

### Why auto-detect CI

The CI-led use case shouldn't require teams to discover and add `--json --summary` to their pipeline definitions. Auto-detection via `CI=true` (universal) and `GITHUB_ACTIONS=true` (GitHub-specific annotations) makes the common case Just Work. The `AGENT=1` / `CLAUDECODE=1` detection extends the same logic to agentic harnesses that benefit from terse machine-parseable output. Explicit flags always override.

## Model Gaps Identified

The ADR depends on the `v1beta1` lineage model. The following items are *not* in the current model and are named here as follow-ups (none block v1 of validate):

1. **`ValueOutput`** — value emits / primitive return channels are not recorded. Validate's reach is currently `FileOutput` + workflow-level `Parameter`s.
2. **Runtime fingerprint completeness** — container digests are sometimes recorded; conda env hashes, executor identity, and base-image provenance are inconsistent.
3. **Workflow `-outputDir` recording** — needed for D14's output join. May already be present in `WorkflowRun.config`; needs confirmation and a stable accessor.
4. **Tagged / named runs** — no first-class alias mechanism (e.g., `lid://golden/v1.0`). Required for the deferred tagged-baseline workflow.
5. **`Workflow.containerImages` / digest pinning** — needs to be stable to support strict runtime mode.

Each gap gets a tracking issue when this ADR is approved.

## Future Work

- **Multi-baseline (N-to-1)** — `--against` repeatable + `--quorum any|all`. Out of v1.
- **`lineage explain <lid>`** — sibling command for deep DAG / causality walks.
- **`lineage tag <lid> <name>`** — alias management for named baselines.
- **`stimulus://` resolver** — a Stimulus-contributed resolver that hydrates fixtures into the validator (validates the SPI's generality).
- **Threshold / tolerance** — bounded checksum sets, file-size deltas. Not in v1; requires careful semantics design.
- **Type-aware comparators** — VCF, BAM, JSON-aware diffs as plugins. Out of scope until D3's checksum-only model proves insufficient.

## Open Questions

These are non-blocking but worth resolving before approval:

1. **Snapshot file extension** — `.json`, `.lineage.json`, `.nflineage`? Discoverable patterns matter for IDE plugins and Git attributes.
2. **`--save-snapshot` default path** — implicit `./<runName>.lineage.json` or require an explicit path?
3. **JUnit / SARIF output formats** — useful for CI test reporters but adds dependencies; defer to first user request?
4. **Confluence with `nextflow lineage diff`** — should there be a non-asserting `diff` sibling, or is `validate --no-fail` the answer?
5. **Permissions** — do we need to assert that the caller can read both lineage records (relevant for shared/network LinStores)?

## Links

- Implementation:
  - CLI: `modules/nextflow/src/main/groovy/nextflow/cli/CmdLineage.groovy` (CmdValidate)
  - Core: `modules/nf-lineage/src/main/nextflow/lineage/cli/LinCommandImpl.groovy#validate`
  - Normalizer: `modules/nf-lineage/src/main/nextflow/lineage/LinNormalizer.groovy`
  - Snapshotter: `modules/nf-lineage/src/main/nextflow/lineage/test/LineageSnapshotter.groovy`
  - Integration tests: `modules/nf-lineage/src/test/nextflow/lineage/LinValidateIntegrationTest.groovy`
- Related: [seqeralabs/stimulus#535](https://github.com/seqeralabs/stimulus/pull/535) (`seqera-compare-runs` skill — Platform + lineage + tasks join)
- Plugin precedent: `adr/20250922-plugin-spec.md`

## More information

- [What is an ADR and why should you use them](https://github.com/thomvaill/log4brains/tree/master#-what-is-an-adr-and-why-should-you-use-them)
- [ADR GitHub organization](https://adr.github.io/)
