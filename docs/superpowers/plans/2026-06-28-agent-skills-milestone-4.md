# Nextflow Agent — Milestone 4 (Agent Skills) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: use superpowers:subagent-driven-development or
> superpowers:executing-plans. Steps use checkbox (`- [ ]`) syntax. Each task is TDD: failing test →
> implement → pass → `git commit -s`. Design: `docs/superpowers/specs/2026-06-28-agent-skills-m4-design.md`.

**Goal:** Add a `skills` agent directive (parallel to `tools`) that lets an agent use Anthropic-style
Agent Skills (`SKILL.md` folders), resolved **locally** from `<baseDir>/skills/<name>/` or **remotely**
from a GitHub repo cached into that same `skills/` dir. Skills are advertised to the model via
langchain4j-skills Tool Mode (`activate_skill` / `read_skill_resource`). Upgrade langchain4j 1.16.3 →
1.17.0. Ship with a working end-to-end POC skill.

**Architecture:** Skills cross the core↔plugin SPI as a portable, langchain4j-free `SkillDescriptor`
(name, description, content, resources) — exactly like `ToolDescriptor`. Core (`SkillResolver`) does all
filesystem/SCM work: ref disambiguation, local `skills/` scan, hand-rolled `SKILL.md` frontmatter parse,
hardened resource load, and remote JGit shallow-clone into a rev-keyed cache dir. The plugin
(`SkillAdapter`) maps descriptors → langchain4j `Skill`s and builds `Skills.from(...)`; the runner attaches
`skills.toolProvider()` additively next to the existing `.tools(map)` and injects
`skills.formatAvailableSkills()` into the system message so the model knows which skills exist.

**Tech Stack:** Groovy 4.0.x (`@CompileStatic`/`@CompileDynamic`), Java (nf-lang), Spock, Gradle, JGit
7.1.1 + SnakeYAML 2.2 (core deps), langchain4j **1.17.0** + `langchain4j-skills:1.17.0-beta27` (plugin).

## Global Constraints
- **Core (`modules/nextflow`, `nf-lang`) MUST NOT import `dev.langchain4j.*`; plugin (`nf-agent`) MUST NOT
  touch `ProcessDef`/`Channel`/`Path`.** Skills cross as portable `SkillDescriptor`/`SkillResource` Maps/DTOs.
- **Never call `executeToolsConcurrently()`.**
- **Tool Mode only** — no `langchain4j-experimental-skills-shell`, no skill-scoped tool/script execution.
- **Remote skill content is untrusted** — resource scanner skips symlinks/`.git` and rejects path escapes;
  caps 64 files / 256 KB; recommend SHA pins. Remote clones are public-repo, single-run for the POC.
- All commits use `git commit -s` (DCO).

---

## File Structure

| File | Status | Responsibility |
|------|--------|----------------|
| `plugins/nf-agent/build.gradle` | Modify | langchain4j → 1.17.0; add `langchain4j-skills:1.17.0-beta27` |
| `modules/nextflow/.../agent/SkillResource.groovy` | **Create** | portable DTO `{relativePath, content}` |
| `modules/nextflow/.../agent/SkillDescriptor.groovy` | **Create** | portable DTO `{name, description, content, resources}` |
| `modules/nextflow/.../agent/AgentRunnerRequest.groovy` | Modify | add `List skills` field |
| `modules/nextflow/.../agent/SkillResolver.groovy` | **Create** | ref disambiguation, local scan, frontmatter parse, resource load, remote clone+cache |
| `modules/nextflow/.../agent/AgentConfig.groovy` | Modify | add `skillsDir` option |
| `modules/nextflow/.../script/AgentBuilder.groovy` | Modify | add `'skills'` to `DIRECTIVES` |
| `modules/nf-lang/.../script/dsl/AgentDsl.java` | Modify | add `void skills(Object... values)` |
| `modules/nextflow/.../script/AgentDef.groovy` | Modify | `getSkills()`; resolve skills (once, via `ownerBaseDir`); set `request.skills`; extend structured guard |
| `plugins/nf-agent/.../agent/SkillAdapter.groovy` | **Create** | `SkillDescriptor`→`Skill`; build `Skills` |
| `plugins/nf-agent/.../agent/LangChainAgentRunner.groovy` | Modify | `isToolRun` gate; attach skills provider; inject catalog; `.tools` only when non-empty |
| `examples/agents/skills/**` | **Create** | POC skill + `main.nf` + README |
| `docs/agent.mdx` | Modify | document `skills` + remote trust note |
| Tests (per task) | Create/Modify | resolver, adapter, request, AgentDef, runner, loading, gated E2E |

---

## Task 0: langchain4j 1.17.0 upgrade + compile-check (gate FIRST)
**Why first:** flush any 1.16→1.17 API breakage in `ChatModelFactory`/`LangChainAgentRunner` before
building skills (the open-ai/core builder methods, esp. `returnThinking`, are unverified from jars).

**Files:** `plugins/nf-agent/build.gradle`.

### Steps
- [ ] Bump in `build.gradle:53-54`: `langchain4j-open-ai:1.17.0`, `langchain4j:1.17.0`; add
  `api 'dev.langchain4j:langchain4j-skills:1.17.0-beta27'`.
- [ ] Compile + run the existing agent gate (no behavior change expected):
  ```bash
  ./gradlew :plugins:nf-agent:compileGroovy
  NXF_SMOKE=1 ./gradlew :plugins:nf-agent:test
  ```
  Expect BUILD SUCCESSFUL. If a builder method moved (e.g. `returnThinking`, `responseFormat`,
  `strictJsonSchema`, the cap-overflow exception type at `LangChainAgentRunner.groovy:182-188`), fix to
  the 1.17.0 signature and note it.
- [ ] Confirm `dev.langchain4j.skills.{Skills,Skill,SkillResource,FileSystemSkillLoader}` resolve on the
  plugin classpath (a throwaway import in a scratch file or `javap`-via-gradle dependency).
- [ ] **Commit:** `git commit -s -m "Bump nf-agent langchain4j to 1.17.0 + add langchain4j-skills"`

---

## Task 1: portable `SkillDescriptor` / `SkillResource` + `AgentRunnerRequest.skills`
**Files:** create `SkillResource.groovy`, `SkillDescriptor.groovy`; modify `AgentRunnerRequest.groovy`;
update the named-args construction test.

### Steps
- [ ] **Failing test:** extend the `AgentRunnerRequest` named-args test (the one asserting all fields) to
  set `skills: [new SkillDescriptor('n','d','c',[])]` and assert it round-trips. Run → FAIL (field missing).
- [ ] Create `SkillResource.groovy` (`@Canonical @CompileStatic`, Apache header): `String relativePath;
  String content`.
- [ ] Create `SkillDescriptor.groovy` (`@Canonical @CompileStatic`, Apache header): `String name; String
  description; String content; List<SkillResource> resources`. Javadoc: portable, langchain4j-free, mirrors
  `ToolDescriptor`.
- [ ] Add `List skills` to `AgentRunnerRequest` (append after `trace`); update the field-order javadoc.
- [ ] Run → PASS; `./gradlew :nextflow:compileGroovy`.
- [ ] **Commit:** `git commit -s -m "feat: portable SkillDescriptor + AgentRunnerRequest.skills"`

---

## Task 2: `SkillResolver` — local resolution, frontmatter parse, resource load
**Files:** create `SkillResolver.groovy`; test `SkillResolverTest.groovy`.

**Interfaces:**
- `static SkillDescriptor SkillResolver.parseSkillDir(Path skillDir)` — parse `SKILL.md` + load resources.
- `static List<SkillDescriptor> SkillResolver.loadLocal(Path baseDir, String name)` — resolve
  `<baseDir>/skills/<name>/` (a single skill dir) or scan it for sub-skills; error if no `SKILL.md`.
- `static Map parseFrontmatter(String text)` — `[name:, description:, content:]` (hand-rolled split).

### Steps
- [ ] **Failing tests** in `SkillResolverTest.groovy` (`@TempDir`):
  - frontmatter: valid `---\nname: x\ndescription: y\n---\nBODY` → `[name:'x',description:'y',content:'BODY']`;
    tolerate leading BOM/blank lines + CRLF; malformed (no frontmatter / unterminated `---` / missing
    name) → throws.
  - resource load: a skill dir with `SKILL.md` + `references/a.txt` → descriptor with one resource
    `relativePath=='references/a.txt'`; a `.git/` dir and a symlink are SKIPPED; >64 files or >256 KB caps.
  - local resolution: `<base>/skills/foo/SKILL.md` → one descriptor named per frontmatter; missing dir →
    `ScriptRuntimeException`.
- [ ] Run → FAIL (class missing).
- [ ] Implement `SkillResolver` (`@CompileStatic`, Apache header), core-only (no langchain4j):
  - `parseFrontmatter`: trim BOM; require first non-blank line `== '---'`; read to next `'---'`;
    `new org.yaml.snakeyaml.Yaml().load(middle) as Map`; `content` = remainder (preserve as-is). Validate
    non-blank `name`/`description`.
  - resource walk: `Files.walk(skillDir)`, filter regular files, exclude `SKILL.md`, exclude any path with
    a `.git` segment, `Files.isSymbolicLink` skip, assert `realPath` stays within `skillDir.realPath`
    (reject escape), enforce 64-file / 256-KB caps (warn + stop). `relativePath` = `skillDir.relativize(p)`.
  - `loadLocal`: `dir = baseDir.resolve('skills').resolve(name)`; if `dir/SKILL.md` exists → one skill; else
    scan immediate subdirs for `SKILL.md`; if none → error.
- [ ] Run → PASS; compile.
- [ ] **Commit:** `git commit -s -m "feat: SkillResolver local resolution + SKILL.md frontmatter parse"`

---

## Task 3: `SkillResolver` — remote GitHub fetch + rev-keyed cache
**Files:** modify `SkillResolver.groovy`; extend `SkillResolverTest.groovy`.

**Interfaces:**
- `static boolean SkillResolver.isRemoteRef(String ref)` — true for `https://github.com/o/r[@rev]`,
  `git@github.com:o/r[@rev]`, `github.com/o/r[@rev]`.
- `static List<SkillDescriptor> SkillResolver.loadRemote(Path baseDir, String ref)` — clone+cache, then
  parse like local.

### Steps
- [ ] **Failing test (OFFLINE, `file://`):** `Git.init()` a temp "remote" repo containing `SKILL.md`
  (+`references/r.txt`), commit it. Call `loadRemote(base, "file://${remote}/.git")` (extend `isRemoteRef`
  test-side or add an overload that accepts a resolved URL) → assert a `SkillDescriptor` is produced and a
  `skills/<name>/` (or `@rev`) cache dir now exists; calling again reuses the cache (no re-clone — assert
  by mtime or a sentinel). Pin `@<sha>` → cache dir suffixed with the rev and HEAD == sha. Model on
  `LocalRepositoryProviderTest.groovy:44,80`.
- [ ] Run → FAIL.
- [ ] Implement: `isRemoteRef` (regex on the three host forms); parse `@rev`; map ref → clone URL + repo
  name; cache dir `baseDir/skills/<repo>[@<rev>]`; if exists → reuse; else
  `Git.cloneRepository().setURI(url).setDirectory(tmp).setDepth(rev?0:1).call().close()`, checkout `rev`
  when set, then atomic-rename tmp → cache dir; scan for `SKILL.md`(root or subdirs) via Task-2 logic;
  zero `SKILL.md` → error. Public/no-credentials for POC.
- [ ] Run → PASS; compile.
- [ ] **Commit:** `git commit -s -m "feat: SkillResolver remote GitHub fetch + rev-keyed local cache"`

---

## Task 4: `skills` directive + `AgentDef` resolution + structured-output guard
**Files:** modify `AgentBuilder.groovy`, `AgentDsl.java`, `AgentConfig.groovy`, `AgentDef.groovy`; tests in
the agent script-loading test + an `AgentDef` resolution test.

### Steps
- [ ] **Failing tests:**
  - loading: an agent with `skills 'foo'` + a local `skills/foo/SKILL.md` LOADS and `AgentDef.getSkills()
    == ['foo']`.
  - resolution: `AgentDef.run()` produces a request whose `skills` has one descriptor named `foo`.
  - guard: an agent with `skills 'foo'` AND a `Record` output type throws `ScriptRuntimeException`
    (skills + structured output unsupported), same as the `tools` case.
- [ ] Run → FAIL.
- [ ] Add `'skills'` to `AgentBuilder.DIRECTIVES`. Add `void skills(Object... values)` (`@Description`) to
  `AgentDsl.DirectiveDsl`. Add `AgentConfig.skillsDir`. Add `AgentDef.getSkills()` (mirror `getTools`).
- [ ] In `AgentDef.run()`: resolve skills ONCE (before the mapper, like `toolSpecs`) via `SkillResolver`
  using `ownerBaseDir(meta, session)` (and `agentConfig.skillsDir` override); set `skills:` on the request
  (independent of the bridge). Extend the guard at line 171 to
  `if( (this.tools || this.skills) && structured ) throw ...` with a message naming skills.
- [ ] Run → PASS; compile core + nf-lang.
- [ ] **Commit:** `git commit -s -m "feat: skills agent directive + AgentDef resolution"`

---

## Task 5: plugin `SkillAdapter` + runner wiring
**Files:** create `SkillAdapter.groovy`; modify `LangChainAgentRunner.groovy`; tests `SkillAdapterTest`,
skills-only runner unit.

### Steps
- [ ] **Failing tests:**
  - `SkillAdapterTest`: `SkillAdapter.toSkills([descriptor])` returns a `Skills` whose
    `formatAvailableSkills()` mentions the skill name/description and whose `toolProvider()` is non-null.
  - runner: a request with `skills` (and empty `toolSpecs`) takes the `runWithTools` path (assert via a
    seam — e.g. that the system message includes the `<available_skills>` block); model can be a stub/mock.
- [ ] Run → FAIL.
- [ ] Create `SkillAdapter` (plugin): map `SkillDescriptor`→`Skill.builder().name().description().content()
  .resources(...)` (resources→`SkillResource.builder()`), return `Skills.from(list)` (keep the handle).
- [ ] In `LangChainAgentRunner`: `run()` gate → `(request.toolSpecs || request.skills) ? runWithTools :
  runSingleShot`. In `runWithTools`: `isToolRun = toolSpecs || skills`; call `.tools(tools)` only when the
  map is non-empty; when `request.skills`, `def skills = SkillAdapter.toSkills(request.skills)`,
  `.toolProvider(skills.toolProvider())`, and append `skills.formatAvailableSkills()` to the composed
  system message. Key the trace-narration nudge off `isToolRun`.
- [ ] Run → PASS; `./gradlew :plugins:nf-agent:compileGroovy`.
- [ ] **Boundary greps:** core has no `dev.langchain4j`; plugin has no `java.nio.file.Path`/`ProcessDef`/`Channel`.
- [ ] **Commit:** `git commit -s -m "feat: nf-agent SkillAdapter + skills toolProvider wiring"`

---

## Task 6: POC skill + example + docs
**Files:** create `examples/agents/skills/{main.nf,README.md,skills/<name>/SKILL.md}`; modify `docs/agent.mdx`.

### Steps
- [ ] Author a small local POC skill whose instructions visibly shape output (e.g. a fixed report format
  or a domain convention the model wouldn't otherwise apply). `main.nf` declares `model`, `instruction`,
  `skills '<name>'`, one input, one `String` output, a prompt that triggers activation (shape from
  `examples/agents/tool/main.nf`). README explains running it (incl. `OPENAI_API_KEY`).
- [ ] Confirm the example LOADS offline (a loading test or `nextflow inspect`).
- [ ] `docs/agent.mdx`: add a `skills` Directives-table row + a Skills section (local `skills/` dir,
  remote `github.com/org/repo[@sha]` form, Tool Mode, **remote content is untrusted — pin a SHA**).
- [ ] **Commit:** `git commit -s -m "docs+examples: agent skills POC + documentation"`

---

## Task 7: Milestone 4 gate
**Files:** none — verification only.

### Exit-criteria steps
- [ ] Boundary: `grep -RIl 'dev\.langchain4j' modules/nextflow/src/main modules/nf-lang/src/main` → empty.
- [ ] Boundary: plugin has no `java.nio.file.Path` / `ProcessDef` / `Channel` imports.
- [ ] `grep -RIn 'executeToolsConcurrently' plugins/nf-agent/src/main` → empty.
- [ ] New unit tests green:
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.SkillResolverTest" --tests "nextflow.agent.AgentRunnerRequestTest" --tests "nextflow.script.*Agent*"
  NXF_SMOKE=1 ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.SkillAdapterTest" --tests "nextflow.agent.*Skill*"
  ```
- [ ] Regression: `./gradlew :nextflow:test --tests "nextflow.agent.*"` and `NXF_SMOKE=1 ./gradlew :plugins:nf-agent:test` → no new failures (the 10 pre-existing M3 SpockTimeout flakes excepted).
- [ ] **Gated E2E (POC):** with `OPENAI_API_KEY`, run the POC agent (or a gated test) end-to-end; confirm
  the trace shows an `activate_skill` call and the output reflects the skill's instructions.

### Milestone 4 exit criteria (all must hold)
1. `skills` is a first-class agent directive; values are local names or `github.com/...` refs.
2. Local skills resolve from `<baseDir>/skills/<name>/`; remote skills clone into a rev-keyed cache under
   the same `skills/` dir and reuse on cache hit.
3. Skills reach the model via langchain4j-skills Tool Mode with the catalog injected into the system
   message; `activate_skill`/`read_skill_resource` work; the POC skill activates end-to-end.
4. Core imports no `dev.langchain4j.*`; plugin touches no `Path`; `executeToolsConcurrently` never called.
5. Skills + structured-output rejected up front; skills-only agents run (null bridge) and terminate cleanly.
6. langchain4j on 1.17.0; new unit tests green; agent suites green; gated E2E green when credentials present.
