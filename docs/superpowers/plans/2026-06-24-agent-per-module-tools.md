# Per-Module Agent Tools — Restore Enforced Schemas (replace aggregate `module_run`)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]`.
>
> **Status: implemented** on branch `agent-aiservices-m3-tools-sandbox` (see commits at the end).

**Goal:** Keep the `tools 'module_run'` directive (the capability that exposes the
script's `include`d modules as agent tools, discovered from `include` statements) but
advertise **each included module as its OWN langchain4j/OpenAI tool** whose function
`parameters` schema IS that module's flattened input schema (required fields,
`additionalProperties:false`, the nf-core `meta.id` convention). This **replaces** the
single aggregate `module_run(module, args)` tool that an earlier draft of Milestone 3
shipped.

## Rationale (the reliability argument)

- A single `module_run` tool must use a generic `args:{type:object,
  additionalProperties:true}`; OpenAI function-calling **cannot enforce field names
  inside a free object**, so the per-module schema was demoted to description prose.
  In practice the model **hallucinated input names** — it sent `reads` though every
  metadata source (the process `path(fastq)`, the `meta.yml`, and the registry) names
  the input `fastq` — and the marshaller silently defaulted the missing required path
  to `[]` → the module ran with empty input and failed.
- The per-module approach puts each module's schema in the tool's function
  `parameters`, so OpenAI **generates/validates** the call against it (required
  `fastq`, `additionalProperties:false`) — the model cannot omit/rename fields. That
  enforcement is the reliability mechanism and is **only achievable with one tool per
  module**.
- **Trade accepted:** more tools advertised (one per included module) vs one
  aggregate. For typical agents (e.g. `isolate-triage` = 3 modules) this is
  negligible.

See the matching *Design note: per-module tools vs single `module_run`* in
[`../specs/2026-06-22-agent-langchain4j-refactor-design.md`](../specs/2026-06-22-agent-langchain4j-refactor-design.md).

---

## Task 1 — Un-suppress the per-module descriptors in `ModuleToolBridge`
**File:** `modules/nextflow/.../agent/ModuleToolBridge.groovy`
- [x] In `wireScalar` and `wireSpec`, register **each** wired module's own
  `ToolDescriptor` (its flattened input schema) in `descriptors` — the per-module
  schema is the tool's function `parameters`. (Previously these descriptors were
  collapsed into one aggregate `module_run` descriptor.)
- [x] Keep the descriptor source precedence: `ModuleMetadataToolSchema` (registry
  metadata, nf-core `meta.id`) when available, else `ModuleSpecToolSchema` (sibling
  `meta.yml`), else the scalar `ProcessToolSchema`.

## Task 2 — Remove the aggregate-tool machinery
**File:** `modules/nextflow/.../agent/ModuleToolBridge.groovy` (+ delete `ModuleRunToolSchema.groovy`)
- [x] Delete the `ModuleRunToolSchema` factory class (no aggregate `module_run`
  descriptor; no `module` enum / generic `args`).
- [x] Remove `callModuleRun(argsJson)` and its `{module, args}` unwrapping. `call()`
  now routes by **tool name**: the `filesystem` tool name → `callFilesystem`; any
  other name → look up the pre-wired module and dispatch via the existing
  `callSpec`/`callScalar` (args ARE the module's flattened inputs, marshalled by
  `ProcessEntryHandler.getProcessArguments`).
- [x] Remove `buildModuleHint` (the per-module summary text that fed the aggregate
  tool's description prose) — the per-module schema now carries that information as
  enforced `parameters`.

## Task 3 — Relocate the module-output read whitelist
**File:** `modules/nextflow/.../agent/ModuleToolBridge.groovy`
- [x] After a non-error module result in `call()`, scan the serialized result for
  absolute path strings and add their parent dirs to the per-record dispatch context's
  readable-dirs whitelist (`whitelistOutputDirs` / `collectPathsFromValue`), so the
  `filesystem` tool may read module outputs. This logic was previously inside
  `callModuleRun`; it now lives on the generic `call()` path keyed off any module tool
  result. Error results (`{"error":…}`) are skipped (no whitelist widening from an
  error message path).

## Task 4 — Keep the metadata feeding per-module schemas
**File:** `modules/nextflow/.../script/AgentDef.groovy`
- [x] `createToolBridge()` still discovers candidate modules from `include` statements
  via `ScriptMeta.getProcessNames()` (all in-scope processes — local + included) when
  `'module_run'` is declared, and still feeds each per-module descriptor from the
  registry `ModuleMetadata` (preferred) or the sibling `meta.yml` `ModuleSpec`
  (fallback). No aggregate descriptor is built; each wired module contributes its own
  `ToolDescriptor`.
- [x] Keep `MODULE_RUN_CAPABILITY = 'module_run'` and the `'filesystem'` capability
  string; keep the per-record sandbox allocation + `DispatchContext` threading.

## Task 5 — Tests + examples
- [x] Bridge tests assert one descriptor **per** wired module (each with its flattened
  input schema), routing by tool name to the right module, marshalling via the existing
  rules, `{"error":…}` on unknown tool name / malformed args.
- [x] `AgentDef` tool-resolution test asserts `tools 'module_run','filesystem'`
  registers the `filesystem` tool plus one descriptor per included module.
- [x] Examples updated to the per-module-tool model (`docs+examples` commit).

---

## Self-review
- The `tools 'module_run'` **directive/capability is unchanged** — only the *shape of
  what it advertises* changes (one tool per module instead of one aggregate).
- Enforcement (the load-bearing property) is restored: each module's flattened input
  schema is the tool's function `parameters`, so OpenAI cannot omit/rename fields — the
  fix for the observed `reads`-vs-`fastq` hallucination.
- Discovery, registry/`meta.yml`/`meta.id` feeding, native dataflow execution, the
  per-record sandbox + output-dir whitelist, and serialized dispatch are all preserved.
- No `ModuleRunToolSchema`, no `callModuleRun`, no `buildModuleHint`, no
  `ToolDescriptor.type` discriminator remain.

## Commits (this branch)
- `feat: module_run advertises full per-module input schema (types, descriptions, nf-core meta.id) to the LLM`
- `feat: module_run fetches registry ModuleMetadata for included modules (prefer registry, fallback meta.yml)`
- `feat: per-module agent tools with enforced schemas (replace module_run aggregate)`
- `fix: per-module tools review findings`
- `docs+examples: per-module agent tools (replace module_run aggregate)`
