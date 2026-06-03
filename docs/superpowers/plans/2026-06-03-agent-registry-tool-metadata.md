# Registry-Sourced Tool Metadata + High-Level Instructions (Plan H)

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Steps use `- [ ]`.

**Goal:** Make a module tool **self-describing to the LLM from the module registry metadata** so the agent `instruction` can be a high-level functional goal. Reuse `module run`'s param→channel binding (`ProcessEntryHandler`) and its nf-core `meta.id` convention, source the tool description + input schema from the registry `ModuleMetadata` (`/api/v1/modules/{name}`, public), and fix the nf-agent `JsonSchemaMapper` so per-field descriptions actually reach the LLM.

**Verified facts (research + adversarial critique, 2026-06-03):**
- `GET /api/v1/modules/{name}` is `@Secured(IS_ANONYMOUS)` (public, no token). `RegistryClient.getModule(name).latest.metadata` → `ModuleMetadata` (npr-client 0.22.0): `getInput():List<ModuleChannel>`, `getOutput():Map<String,ModuleChannel>` (unordered → NOT usable for positional marshalling), `getTools():List<ModuleTool>` (name/description/homepage/documentation URIs), `getDescription()`. `ModuleChannelItem`: `getName/getType/getDescription/getPattern/get_enum` (getters can return null — must null-guard).
- `module run` binds params via `ProcessEntryHandler.getProcessArguments(processDef, params)` (`ProcessEntryHandler.groovy:156`): dot-params → nested maps; per-input value lookup by name; type coercion from the sibling `meta.yml` (`file`→`Nextflow.file`, `map`→Map, etc.); tuple assembly → the channel value `[[id:..], file(..)]`. The `meta.id` is a **hardcoded nf-core convention** in `CmdModuleView.inferNfCoreParam` (`CmdModuleView.groovy:233`), NOT a sub-schema.
- **SHOWSTOPPER:** `JsonSchemaMapper.toElement`/`toObjectSchema` (`plugins/nf-agent/.../JsonSchemaMapper.groovy`) switch on `type` only and **drop per-property `description` and ignore `enum`** — so the rich per-field hints never reach the LLM through `inputSchema` today. Only the tool-level `ToolSpecification.description` survives (and that path works).
- Marshalling must stay on ordered structure (ProcessDef + ordered `ModuleSpec`); the registry `output` Map is unordered. Descriptor (from `ModuleMetadata`) and execution (from ProcessDef/meta.yml) are two views of the same `meta.yml` → add a consistency check.

---

## Task 1 — Fix `JsonSchemaMapper` to propagate per-property descriptions (+ enum) — highest leverage, independent
**Files:** `plugins/nf-agent/src/main/nextflow/agent/JsonSchemaMapper.groovy`; test.
- [ ] Failing test (`JsonSchemaMapperTest`): a portable schema with a property carrying a `description` (and a `meta` object with nested `properties:{id:{type:string}}`) → the built `JsonObjectSchema` exposes that property `description` and the nested `id` property. Also an `enum` fragment → `JsonEnumSchema`/enum carried.
- [ ] Implement: in `toElement`/`toObjectSchema`, call `.description(fragment.description)` on each langchain4j schema element when present; recurse nested `object` `properties` + `required`; map a `string` with `enum` → `JsonEnumSchema` (or `JsonStringSchema` with enum if the API requires). Verify the langchain4j 1.x builders (`JsonStringSchema`/`JsonObjectSchema`/`JsonNumberSchema`/`JsonEnumSchema`) expose `.description(...)` via `javap`.
- [ ] Run `:plugins:nf-agent:test` → green. Commit `feat: propagate per-field descriptions (and enum) into langchain4j tool schemas`.

## Task 2 — Reuse `ProcessEntryHandler` binding for tool marshalling
**Files:** `modules/nextflow/src/main/groovy/nextflow/script/ProcessEntryHandler.groovy` (make the binding reusable); `nextflow/agent/ModuleToolBridge.groovy`; tests.
- [ ] Factor the param→args binding into a reusable, accessible method: `static List<Object> getProcessArguments(ProcessDef proc, Map params, Path moduleSpecPath=null)` (or extract a `ProcessInputBinder`), preserving the existing `module run` behavior (dot-params, type coercion from meta.yml, tuple assembly, `Nextflow.file` for files). Keep `module run`'s existing call working unchanged.
- [ ] Rewire `ModuleToolBridge` (spec-driven path) input marshalling: at dispatch, parse the LLM args JSON → a params-style Map; `args = ProcessEntryHandler.getProcessArguments(toolProcessDef, argsMap, specPath)`; bind `args[i]` onto the pre-wired `toolIns[i]`. Replaces `marshalInput`. Output serialization (ordered `ModuleSpec`, eval-skip) is unchanged. The number of input queues = the process's declared input-channel count.
- [ ] Test: a tool dispatch with `{meta:{id:'s1'}, reads:'/abs/x.fq'}` (or flattened `meta.id`) builds the same channel value `[[id:'s1'], file(...)]` and the real process runs (reuse the `AgentModuleSpecToolTest` harness). Confirm `module run` tests still pass.
- [ ] Commit `feat: reuse ProcessEntryHandler param binding for agent tool marshalling`.

## Task 3 — Source tool description + input schema from registry `ModuleMetadata` (with `meta.id` convention)
**Files:** `AgentDef.groovy` (`resolveRegistryTool`/`createToolBridge`); new `nextflow/agent/ModuleMetadataToolSchema.groovy`; `ModuleToolBridge.groovy` (ctor + `wireSpec`); tests.
- [ ] `resolveRegistryTool`: reuse the `RegistryClient` already built to `getModule(ref.fullName).latest?.metadata` (try/catch → `log.warn` → null on failure). Thread the `ModuleMetadata` to the bridge alongside the ProcessDef + the sibling `ModuleSpec`. (Use the install-resolved version if cheaply available; else `.latest` + the consistency check below.)
- [ ] `ModuleMetadataToolSchema`: `inputSchema(ModuleMetadata)` mirroring `ModuleSpecToolSchema` but fed from `ModuleChannelItem` (null-guard `getName`), folding in `description`, `pattern` (append to description), `enum`. For an nf-core `meta` map item, emit a nested object `{type:object, description:<hint>, properties:{id:{type:string,description:'sample identifier'}}, additionalProperties:true}` — mirroring `inferNfCoreParam`'s `meta.id` convention. `buildDescription(ModuleMetadata)` = module description + `tools[]` (name/description/homepage/documentation) + output-shape prose (emit names preserved).
- [ ] `ModuleToolBridge.wireSpec`: when `metadata != null`, build `ToolDescriptor.description`/`inputSchema` from `ModuleMetadataToolSchema`; else fall back to `ModuleSpecToolSchema` (offline/local-file/in-scope). Marshalling/output unchanged (ProcessDef + `ModuleSpec`).
- [ ] Consistency check: at wire time, warn if the flattened input property names from `ModuleMetadata` ≠ from `ModuleSpec` (silent-drift guard).
- [ ] Tests: a registry tool (use the `AgentRegistryToolTest` local-install harness; add a `meta.yml` + a stand-in `ModuleMetadata` if the install path can supply one, or unit-test `ModuleMetadataToolSchema` directly) → descriptor description contains the module/tool prose; inputSchema has `meta.properties.id` + per-field descriptions. Confirm the offline fallback still works.
- [ ] Commit `feat: source agent tool description + schema from registry module metadata`.

## Task 4 — High-level instructions + verification + real run
- [ ] Simplify `examples/agents/skesa/main.nf` + `isolate-triage/main.nf` instructions to functional goals (drop "call nf-core/skesa", "build meta as a map with id", etc.) — rely on the registry-sourced tool description/schema.
- [ ] Full keyless agent suite green (core + nf-lang + plugin); gated E2E skipped.
- [ ] Real `./launch.sh` run of the skesa example with the **high-level** instruction → confirm the LLM still calls skesa with `meta={id:..}` (from the schema, not the instruction) and `ASSEMBLY=...`. Update the ADR.
- [ ] Commit `docs: high-level agent instructions; mark registry-sourced tool metadata delivered`.

## Open/deferred
- Version-correct metadata fetch (surface the resolved version from `ModuleResolver.resolve`) — use `.latest` + consistency-warning for v1; note as follow-up.
- In-session metadata cache — none for v1 (one GET per tool at wiring).

## Self-review
- Reuses `module run`'s `ProcessEntryHandler` binding + `meta.id` convention (no parallel convention invented). Descriptor from registry `ModuleMetadata`; execution from ProcessDef + ordered `ModuleSpec`. The `JsonSchemaMapper` fix is the load-bearing change that makes per-field hints (incl. the `meta` nested schema) reach the LLM. Offline `meta.yml` fallback preserved. Consistency check guards descriptor/executor drift.
