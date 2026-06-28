# Nextflow Agent — Milestone 4: Agent Skills (design)

**Status:** design (reviewed by 3-lens adversarial panel, 2026-06-28) • **Branch:** `agent-aiservices-m4-skills` (stacked on M3 `agent-aiservices-m3-tools-sandbox`)

## Goal

Add a `skills` agent directive (parallel to `tools`) so a Nextflow `agent` can use Anthropic-style
**Agent Skills** — `SKILL.md` folders that progressively disclose instructions (and optional reference
files) to the model. Skills can be **local** (a `skills/<name>/` dir beside the running `main.nf`) or
**remote** (fetched from a GitHub repo and cached into that same local `skills/` dir). Skills extend the
agent's capability like tools do. Ship as milestone M4 with a working end-to-end POC skill. Upgrade
langchain4j 1.16.3 → 1.17.0.

## Background (verified)

langchain4j ships `langchain4j-skills` (Tool Mode). It exposes two LLM tools — `activate_skill`
(returns a skill's full instructions on demand) and `read_skill_resource` (returns a bundled file's
content) — and the model is told which skills exist via a catalog string. **Tool Mode has no filesystem
access at inference and runs no arbitrary code** — only the registered tools execute. (A separate
`langchain4j-experimental-skills-shell` adds an unsandboxed `run_shell_command`; it is **out of scope**
for M4, same rationale as the deferred `bash` tool.)

Verified against the repo and the actual jars (`/tmp/lc1170.jar`, `/tmp/lcskills.jar`):
- **The `.tools` vs `.toolProvider` mutual-exclusion is gone** (already in 1.16.3, confirmed in 1.17.0):
  `ToolService` holds a `Set<ToolProvider>` and separate static-tool fields; `.tools(map)` + any
  providers coexist and merge. The skills `ToolProvider` is **additive** next to the module/filesystem
  `.tools(map)`. Bytecode-confirmed: with a non-empty provider set, `createContextFromStaticToolsAndProviders`
  yields a non-empty tool context even when the static map is empty — so a **skills-only** agent runs.
- `Skills.from(Collection<Skill>) : Skills`; `Skills.toolProvider() : ToolProvider`;
  `Skills.formatAvailableSkills() : String` (an `<available_skills>…</available_skills>` catalog block).
- `Skill.builder().name(s).description(s).content(s).resources(Collection<SkillResource>).build()` and
  `SkillResource.builder().relativePath(s).content(s).build()` — build skills **programmatically** from
  in-memory data; no `Path`, no `FileSystemSkillLoader` in the plugin.
- For Tool-Mode skills with no skill-scoped tools, the skills provider's `isDynamic()` is **false** →
  `provideTools` runs once per service; activation state is derived from chat-message history at
  execute-time (not stored in the provider). Building the provider per record over a fresh `ChatMemory`
  is therefore cheap and leak-free.
- `skills` directive needs **no grammar change**: directive names are gated only by
  `AgentBuilder.DIRECTIVES` (runtime) and `AgentDsl.DirectiveDsl` (resolution scope).
- **JGit** `org.eclipse.jgit:7.1.1` and **SnakeYAML** `2.2` are direct `api` deps of `modules/nextflow`
  (build.gradle:63,66); `Git.cloneRepository()…setDepth(1).call()` is already used at
  `LegacyRepositoryStrategy.groovy:79-89`.

> **Unverifiable from the supplied jars → mandatory compile-check (first impl task):** the
> `OpenAiChatModel.builder()` methods (`apiKey/modelName/timeout/listeners/returnThinking/responseFormat/strictJsonSchema/build`),
> `ResponseFormat`/`ResponseFormatType.JSON`/`JsonSchema`/`ChatModelListener` (langchain4j-core), and the
> exact exception thrown on `maxSequentialToolsInvocations` overflow live in `langchain4j-core`/`-open-ai`
> (not in the jars reviewed). They are **expected** present in 1.17.0 but must be flushed by recompiling
> the plugin + running the nf-agent gate before building anything else. `returnThinking` is the
> highest-churn risk.

## Hard invariants (carried from M1–M3)

- Core (`modules/nextflow`, `nf-lang`) **must not import `dev.langchain4j.*`**.
- The `nf-agent` plugin **must not touch `ProcessDef`/`Channel`/`Path`**.
- Tool dispatch stays sequential on the calling thread; never `executeToolsConcurrently()`.

## Core design decisions

### 1. The `skills` directive
- Add `'skills'` to `AgentBuilder.DIRECTIVES`; add `void skills(Object... values)` (with `@Description`)
  to `AgentDsl.DirectiveDsl`, mirroring `tools`. Add `AgentDef.getSkills()` mirroring `getTools()`.
- Value shape identical to `tools`: `skills 'a'` → `['a']`; `skills 'a','b'` → `['a','b']`.

### 2. Local vs remote reference disambiguation (ordered predicate list)
A skill reference is resolved by this **ordered** list (no overloading of the bare `org/repo` form, which
the registry-module rule already claims — avoids the SEV-1 collision):
1. **Remote** iff it has an explicit GitHub host: `https://github.com/<org>/<repo>[@rev]`,
   `git@github.com:<org>/<repo>[@rev]`, or `github.com/<org>/<repo>[@rev]`. → fetch + cache (decision 4).
2. **Local** otherwise: a name resolving to `<baseDir>/skills/<name>/` (a dir containing `SKILL.md`).
3. Anything else → `ScriptRuntimeException` ("skill `x` not found under skills/ and is not a GitHub ref").
- `<baseDir>` is obtained from the **existing** `ownerBaseDir(meta, session)` helper
  (`AgentDef.groovy:708-715`: `getModuleDir()` → `session.baseDir` → cwd), NOT a fresh nullable
  `getModuleDir()` call. When the agent is itself loaded from a registry/remote module, `<baseDir>` is
  that module's dir — document that `skills/` then resolves there.

### 3. SPI shape — portable `SkillDescriptor` (keeps the invariant)
Skills cross the core→plugin boundary as a portable, langchain4j-free DTO, exactly like `ToolDescriptor`:
- **New core types** `nextflow.agent.SkillDescriptor { String name; String description; String content;
  List<SkillResource> resources }` and `nextflow.agent.SkillResource { String relativePath; String
  content }` (both `@Canonical`, langchain4j-free).
- **New field** `AgentRunnerRequest.skills : List<SkillDescriptor>` (append to the `@Canonical` field
  list — AgentDef constructs with named args so order is irrelevant; the named-args-construction test
  must add the field).
- **Core resolves** each directive entry to one or more skill directories (local, or remote-cloned),
  parses each `SKILL.md`, and eagerly loads bundled resources into `SkillDescriptor`s.
- **`SKILL.md` parsing (hand-rolled frontmatter split — NOT `fromYaml`):** snakeyaml does not split
  frontmatter. Detect a leading `---` line, find the next `---`, `new Yaml().load(...)` ONLY that middle
  block → `name`+`description`; the remainder is `content`. Tolerate a leading BOM / blank lines / CRLF
  so a skill that loads in langchain4j's own `FileSystemSkillLoader` also loads here. Malformed (no
  frontmatter, unterminated `---`, missing `name`/`description`) → clear `ScriptRuntimeException`.
- **Resource loading (hardened):** resources = regular files under the skill dir other than `SKILL.md`,
  scanned with `Files.walk`, **skipping `.git/`, skipping symlinks, rejecting any path escaping the skill
  dir**, capped at **64 files / 256 KB total** (over-cap → skip remainder + warn). Each → `SkillResource`.
- **Plugin maps** each `SkillDescriptor` → `Skill.builder()…build()` (resources → `SkillResource.builder()…`).
  Mirrors `ToolDescriptor`→`ModuleToolAdapter`; no `Path` crosses; core stays langchain4j-free.

### 4. Remote fetch + cache
- Shallow clone via JGit into a **rev-keyed** cache dir (fixes the SEV-1 "pinned rev ignored" bug):
  `<baseDir>/skills/<repo>/` for the default branch, `<baseDir>/skills/<repo>@<rev>/` when `@rev` is
  given; checkout `<rev>` after clone; **close the returned `Git`** (the existing code does).
- **Cache:** if the (rev-keyed) dir exists, reuse as-is (offline-friendly). Clone first into a temp dir
  then atomically rename into place, so a failed/interrupted clone never leaves a half-populated cache.
- **Repo→skill(s):** after fetch, scan for a root `SKILL.md` (whole repo = one skill) or immediate-subdir
  `SKILL.md`s (each a skill) — same scan as a local `skills/` dir. **Zero `SKILL.md` → error.**
- **Duplicate skill `name`s** across all resolved descriptors → error (langchain4j requires unique names).
- Scope notes (follow-ups, not POC-blocking): credentialed/private repos via `nextflow.scm.ProviderConfig`;
  a concurrent-run clone lock; auto-`.gitignore` of `skills/`. POC targets public repos, single run.

### 5. Plugin wiring & tool-loop entry
- `LangChainAgentRunner.run()`: gate becomes `(request.toolSpecs || request.skills) ? runWithTools :
  runSingleShot`. Introduce `boolean isToolRun = request.toolSpecs || request.skills` and key the
  trace-narration nudge (`composeSystemMessage`, currently gated on `toolSpecs`) off `isToolRun`, so a
  skills-only `-with-agent-trace` run still narrates.
- In `runWithTools`: build the native `tools` map as today; **only call `.tools(tools)` when non-empty**;
  when `request.skills` is non-empty, `Skills skills = SkillAdapter.toSkills(request.skills)` →
  `.toolProvider(skills.toolProvider())`.
- **Inject the catalog (BLOCKER fix):** append `skills.formatAvailableSkills()` to the composed
  **system message** so the model knows which skills exist (its `activate_skill` `skill_name` parameter
  is free-text, not an enum — without the catalog the model cannot name a skill and activation never
  fires). Keep the `Skills` handle to call both `toolProvider()` and `formatAvailableSkills()`.
- **Skills-only / null-bridge:** skill descriptors are resolved in `AgentDef.run()` **independently of
  `createToolBridge`** and set on `request.skills` unconditionally; a skills-only agent intentionally has
  `bridge == null` / `toolSpecs == null` and terminates via the normal map-operator path (skills add no
  dataflow nodes, so there is nothing to poison). `applyAgentOperator`'s `bridge?.close()` is already
  null-safe.

### 6. Config additions
`agent.skillsDir` (optional) overrides the default `<baseDir>/skills` location — added now because the
default writes remote clones into the user's source tree. Default unchanged when unset.

### 7. Security posture
Tool Mode only; no shell-mode skills; no skill-scoped tool/script execution wired in M4. `activate_skill`
returns text; `read_skill_resource` returns core-preloaded content — no inference-time filesystem access.
**Remote skill content is untrusted model input** (a remote `SKILL.md` becomes model instructions — a
prompt-injection/instruction-smuggling surface, worsened by moving refs): recommend pinning to a commit
SHA; the resource scanner rejects symlinks/path-escapes and skips `.git/` (decision 3). The clone runs no
repo code. Document this in `docs/agent.mdx`.

### 8. langchain4j 1.17.0 upgrade
- `plugins/nf-agent/build.gradle`: bump `langchain4j` + `langchain4j-open-ai` `1.16.3 → 1.17.0` (stable);
  add `dev.langchain4j:langchain4j-skills:1.17.0-beta27`.
- **First implementation task** = the compile-check above (recompile plugin + run nf-agent gate) to flush
  any 1.16→1.17 breakage in `ChatModelFactory`/`LangChainAgentRunner` BEFORE building skills.

### 9. POC scope & verification
- **Local** POC skill `examples/agents/skills/skills/<name>/SKILL.md` whose instructions visibly change
  the agent's output (so activation is observable), + example `examples/agents/skills/main.nf` declaring
  `skills '<name>'` (shape from `examples/agents/tool/main.nf`) + a README.
- **Unit/offline tests:** frontmatter parse (good + malformed); resource scan (caps, symlink/`.git`
  skip); local skill resolution; **remote fetch via a `file://` local clone** (`Git.init()` a temp repo
  with a `SKILL.md`, clone shallow into a temp `skills/`, assert the `SkillDescriptor` + cache-hit reuse)
  — deterministic, offline, matches `LocalRepositoryProviderTest.groovy:44,80`; NOT a real public-repo
  clone (which is skipped under `NXF_SMOKE` and would make "remote is tested" false). Plugin `SkillAdapter`
  unit (descriptor→Skill→toolProvider). Skills-only runner path unit.
- **Gated E2E** (`@Requires({ System.getenv('OPENAI_API_KEY') })`, per `AgentToolEndToEndTest.groovy:33`):
  run the local POC, assert the trace shows an `activate_skill` call and the answer reflects the skill.

## File-level change set

| File | Action | Responsibility |
|------|--------|----------------|
| `modules/nextflow/.../script/AgentBuilder.groovy` | modify | add `'skills'` to `DIRECTIVES` |
| `modules/nf-lang/.../script/dsl/AgentDsl.java` | modify | add `void skills(Object... values)` (`@Description`) |
| `modules/nextflow/.../agent/SkillDescriptor.groovy` | **create** | portable DTO (name, description, content, resources) |
| `modules/nextflow/.../agent/SkillResource.groovy` | **create** | portable resource DTO (relativePath, content) |
| `modules/nextflow/.../agent/SkillResolver.groovy` | **create** | resolve refs → SkillDescriptors: ref disambiguation, local scan, hand-rolled SKILL.md frontmatter parse, hardened resource load, remote JGit clone+cache (rev-keyed, atomic), dup-name/zero-skill errors |
| `modules/nextflow/.../agent/AgentConfig.groovy` | modify | add `skillsDir` option |
| `modules/nextflow/.../agent/AgentRunnerRequest.groovy` | modify | add `List skills` field |
| `modules/nextflow/.../script/AgentDef.groovy` | modify | `getSkills()`; resolve skills (once, pre-ignition, via `ownerBaseDir`); set `request.skills`; extend structured-output guard to `(tools \|\| skills)` |
| `plugins/nf-agent/build.gradle` | modify | bump 1.17.0 + add langchain4j-skills:1.17.0-beta27 |
| `plugins/nf-agent/.../agent/SkillAdapter.groovy` | **create** | `SkillDescriptor`→`Skill`; `Skills.from(...)` (return the `Skills` handle) |
| `plugins/nf-agent/.../agent/LangChainAgentRunner.groovy` | modify | `isToolRun` gate; attach skills `toolProvider`; inject `formatAvailableSkills()` into system msg; `.tools` only when non-empty |
| `examples/agents/skills/**` | **create** | POC skill + `main.nf` + README |
| `docs/agent.mdx` | modify | document `skills` directive + remote-skill trust/prompt-injection note |
| tests (per task) | create/modify | SkillResolver (parse/scan/local/remote-file://), SkillAdapter, request named-args, AgentDef skills+guard, skills-only runner, loading test, gated E2E |

## Risks / open questions (residual, documented)
- **Beta API (`langchain4j-skills:1.17.0-beta27`) churn** — pinned exact; the compile-check task fails fast if a builder signature differs.
- **Remote content trust** — see §7; SHA-pinning recommended, moving refs reused-on-cache-hit by design.
- **Source-tree writes / concurrency** — remote clones land in `<baseDir>/skills/`; `agent.skillsDir`
  override added; concurrent-run clone lock + auto-`.gitignore` are follow-ups (POC = single run).
- **Resource size caps** (64 files/256 KB) are heuristic; revisit if real skills need more.
