# Inline small structured tool outputs to the LLM (Plan I)

> REQUIRED SUB-SKILL: superpowers:subagent-driven-development.

**Goal:** Let a module tool's **small, structured file output** (e.g. `nf-core/assemblyscan`'s stats JSON) be **read and returned to the LLM as content** so it can reason over the numbers — while bulk/binary data (contigs, BAMs) stays an opaque path handle for chaining. This removes the hand-written scalar QC process and lets real nf-core stats modules drive LLM decisions.

**Decisions (locked):** default cap **32 KB** (config `agent { maxToolOutputInlineSize = '32 KB' }`); JSON returned as **raw text**; **inference-only** (format + size, no per-output annotation); too-big readable output → **path handle + note** (no truncation of structured data). Per-loop cumulative budget = deferred.

## The decision (per file output, at serialization time)
```
ext = file extension (lowercase)
TEXT_LIKE = {json, tsv, csv, txt, tab, yaml, yml, log, md}
if ext ∉ TEXT_LIKE                         → opaque handle (abs path string)   // data/binary, chainable
else if Files.size(file) > maxInlineBytes  → { path: <abs>, note: "content not inlined: <size> exceeds <cap>" }
else if looksBinary(first 8 KB has NUL)    → opaque handle (safety net)
else                                       → file content as a UTF-8 String    // the LLM reads it
```

## Tasks

### Task 1 — the inline capability (core)
- New `nextflow/agent/ToolOutputReader.groovy`: `static Object readOrHandle(Path file, long maxBytes)` implementing the rule above (helpers: `extensionOf`, `looksBinary`, human-readable size). Unit test `ToolOutputReaderTest`: small `.json` → its content; oversized `.json` → `[path, note]`; `.fa` → abs path; a NUL-containing `.txt` → abs path; missing-ext → handle.
- `AgentConfig`: add `@ConfigOption MemoryUnit maxToolOutputInlineSize` (default `32 KB`) + a `long maxToolOutputInlineBytes()` accessor (default 32768). Test in `AgentConfigTest`.
- `AgentDef.run`: read `agentConfig().maxToolOutputInlineBytes()`, pass to the bridge.
- `ModuleToolBridge`: store `maxInlineBytes`; in the FILE branch of `serializeValue` (currently `pathString(path)`), call `ToolOutputReader.readOrHandle(path, maxInlineBytes)` instead. Binary/data and non-file values unchanged. Add an integration assertion (reuse the real-executor harness): a module whose output is a small `.json`/`.txt` file → the dispatch result carries the FILE CONTENTS (not a path); a `.fa` output → still a path. Commit `feat: inline small structured tool outputs to the LLM (size-capped, format-inferred)`.

### Task 2 — rebuild the real-world example around real nf-core/assemblyscan + fix compile + real run
- Replace the hand-written `assembly_qc` Groovy process in `examples/agents/isolate-triage/main.nf` with the registry tool `nf-core/assemblyscan` (its `report` JSON is small + `.json` → inlined → the LLM reads N50/#contigs). This also REMOVES the v2-illegal Groovy (`final`/C-for/`Channel.of`) that broke compilation. Tools become `nf-core/skesa, nf-core/assemblyscan, nf-core/prokka` (all registry, all self-describing). Keep the high-level instruction + QC gate (now reasoning over assemblyscan's inlined JSON). Use `channel.of` (lowercase). Ensure the script compiles under `nextflow.enable.types = true`.
- Verify: full keyless agent suite (in isolation). Then a REAL run (live key + Wave/Docker) demonstrating: skesa assembles → assemblyscan runs → its **JSON stats are inlined to the LLM** → the LLM gates on real N50/#contigs → PASS(annotate)/FAIL verdict. Capture the verbatim tool-call sequence + the inlined-stats decision + the `TRIAGE:` verdict. Update the ADR. Commit.

## Open/deferred
Per-loop cumulative inline budget; format-aware previews (head rows / head-tail) for big tabular/log; `{content, path}` dual form if a readable output ever needs chaining; explicit per-output inline/handle override.

---

## DELIVERED (2026-06-04)

**Task 1 — inline tool outputs:** DONE. `ToolOutputReader` (format+size+binary-sniff rule), `agent.maxToolOutputInlineSize` (32 KB default), bridge wiring; 10 unit tests + integration test (real bridge+executor inlines a small `.json`). Spec + code-quality reviews passed. Commits `b57361063`, `12dd5ee87`, `fb4020a65`, `7f4bbbd0a`, `9f77f5478` (the last keeps `AgentRegistryToolTest` on an opaque-path `.dat` after inlining changed `.txt` behavior).

**Task 2 — example rebuilt + live run:** DONE. `isolate-triage` rebuilt around real `nf-core/{skesa,assemblyscan,prokka}` (removed the buggy Groovy QC process). The live end-to-end run surfaced **two pre-existing agent-feature bugs** (only exposed by running multiple/real multi-output nf-core tools), both fixed:

- **Bug A — versions/eval `topic` output blocks the dispatch** (commit `3c02dde46`). The registry meta.yml types the eval/version component inconsistently (skesa `eval`, assemblyscan `string`), so `isEvalOutput` missed assemblyscan's `versions_assemblyscan` and the bridge blocked forever on its topic-source channel's `.val` → the run hung after the tool call. Fix: skip outputs whose **ProcessDef** `OutParam.getChannelTopicName() != null` (authoritative; meta.yml type is unreliable). Test: `AgentModuleSpecToolTest` "topic-routed versions output typed as string".

- **Bug B — multi-tool container bleed** (commit `f1f44ce8e`). Every bridge-compiled tool module shared `session.binding`, so `BaseScript.setup()`'s `moduleDir` write was last-wins (prokka), making `conda "${moduleDir}/environment.yml"` resolve to one shared container for all tools (skesa ran in prokka's container). Fix: give each tool module its own `ScriptBinding` (seeded with owner params) via a new opt-in `ScriptLoaderV2.setMainBinding`, mirroring `IncludeDef.loadModuleV2`. Test: `AgentModuleDirIsolationTest` (3 tools, distinct `moduleDir`). Root cause verified by a 4-investigator + adversarial workflow (high confidence).

**Live end-to-end proof:** `tools 'nf-core/skesa','nf-core/assemblyscan','nf-core/prokka'`, high-level instruction, gpt-5-mini + Wave/Docker. The LLM called skesa → chained its contigs into assemblyscan → assemblyscan's `isolate_001.tsv` stats were **inlined** → the LLM read real stats (`n50_contig_length=863, total_contig=24`), applied the QC gate (N50 863 < 20000 → FAIL), skipped prokka, and emitted: `TRIAGE: FAIL isolate_001: fragmented (N50=863, contigs=24), needs manual review`. Clean exit (`[SUCCESS] completed=2 failed=0`).

**Known test-infra flakiness (not a product bug):** the FULL `AgentModuleSpecToolTest` class flakes at session shutdown when multiple `topic`-creating agent sessions run back-to-back in one JVM (GPars cross-session topic-operator issue). Each method passes in isolation; the single-session live run terminates cleanly. Run agent test classes individually.
