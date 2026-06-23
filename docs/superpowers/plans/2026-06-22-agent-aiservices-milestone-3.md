# Nextflow Agent Refactor — Milestone 3 (filesystem tool, work-dir sandbox, include-driven module_run) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give each agent invocation its own work-dir sandbox; add a generic `filesystem` tool (read/write/list/exists) confined to it; and expose a single generic `module_run(module, args)` tool whose available modules are inferred from the script's `include` statements — with `tools` becoming a capability list (`'filesystem'`, `'module_run'`).

**Architecture:** `tools` entries are capability strings, validated at agent-run time in `AgentDef.createToolBridge`. Module tools are discovered from `ScriptMeta` imports (the `include`d processes), pre-wired into the dataflow graph as today, and surfaced to the LLM as ONE `module_run` tool (a `module` enum + a generic `args` object); the description aggregates each module's input hints. Generic-tool execution runs inside a per-invocation work-dir sandbox enforced by a new `SandboxGuard`. The per-record context (work dir + accumulated module-output read-whitelist) is threaded core-side via a `ThreadLocal` set by the agent operator around `runner.run(...)` and read by `ModuleToolBridge.call(...)` on the same thread — leveraging M1's invariant that AiServices runs tools sequentially on the calling thread. The `nf-agent` plugin is UNCHANGED.

**Tech Stack:** Groovy 4.0.x (`@CompileStatic`/`@CompileDynamic`), Java (nf-lang), Spock, Gradle (`./gradlew`), langchain4j 1.16.3 (from M1).

## Global Constraints

- **Core `modules/nextflow` must never import `dev.langchain4j.*`; the `nf-agent` plugin must never touch `ProcessDef`/`Channel`/`Path`.** Because the plugin cannot touch `Path`, the per-invocation work dir/whitelist is threaded **core-side via a `ThreadLocal`**, never through the plugin or `AgentRunnerRequest`.
- **Never call `executeToolsConcurrently()`** — tool calls stay sequential on the calling thread (the `ThreadLocal` context relies on this invariant; document it where the context is read).
- The per-record context must **not** live as a mutable field on the shared, pre-ignition `ModuleToolBridge` — it is `ThreadLocal` (thread-scoped, cleared in a `finally`), keeping the bridge stateless across records.
- `module_run` is **one** LLM tool: `{module: string enum of included module names, args: object}`. The internal per-module `Tool`s are keyed by module name and are NOT advertised individually.
- The `filesystem` tool: **writes** confined to the agent work dir; **reads** allowed within the work dir + the per-invocation module-output whitelist; `..`/absolute-outside/symlink-escape rejected by `SandboxGuard`. `file://` work dirs only — return `{"error":...}` otherwise.
- **`bash` and `nextflow run` are deferred** — not in this milestone.
- **Scope decision (documented):** the legacy module-reference string forms in `tools` (`'nf-core/...'`, local paths, in-scope process names) remain working-but-deprecated in M3 to avoid mass test/example breakage; examples are migrated to the new model; full removal of the legacy path is a follow-up. (See "Open question" at the end.)
- Tool dispatch returns the existing `{...}` / `{"error":...}` JSON contract; `call()` stays `synchronized`.
- All commits use `git commit -s` (DCO).

---

## File Structure

| File | Status | Responsibility |
|------|--------|----------------|
| `modules/nextflow/src/main/groovy/nextflow/script/ScriptMeta.groovy` | Modify | Add `getIncludedProcessNames()` to enumerate `include`d process names. |
| `modules/nextflow/src/main/groovy/nextflow/agent/SandboxGuard.groovy` | **Create** | Pure path-containment logic for the filesystem tool. |
| `modules/nextflow/src/main/groovy/nextflow/agent/DispatchContext.groovy` | **Create** | Per-invocation context: work dir + mutable readable-dirs set. |
| `modules/nextflow/src/main/groovy/nextflow/agent/FilesystemToolSchema.groovy` | **Create** | `ToolDescriptor` factory for the `filesystem` tool. |
| `modules/nextflow/src/main/groovy/nextflow/agent/ModuleRunToolSchema.groovy` | **Create** | Build the single `module_run` `ToolDescriptor` (enum + args + aggregated description). |
| `modules/nextflow/src/main/groovy/nextflow/agent/ModuleToolBridge.groovy` | Modify | `ThreadLocal<DispatchContext>`; route `call()` for `filesystem`/`module_run`; `callFilesystem`; `callModuleRun` (delegates to existing `callScalar`/`callSpec` + records output dirs); capability-aware descriptors. |
| `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` | Modify | `createToolBridge`: capability handling + include discovery; `run()`: per-invocation work dir + set/clear the `ThreadLocal`; cleanup on stop. |
| `examples/agents/*` | Modify | Migrate to `include {...}` + `tools 'module_run'`; add a `filesystem` example. |
| `docs/agent.mdx` | Modify | Document the capability model, `module_run`, `filesystem`, include-driven discovery. |
| Tests (per task) | Create/Modify | SandboxGuard, ScriptMeta, bridge generic dispatch, AgentDef resolution, gate. |

---

## Task 1: `ScriptMeta.getIncludedProcessNames()`

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/script/ScriptMeta.groovy` (near `getProcessNames()`/`getLocalProcessNames()`, ~lines 272-295)
- Test: `modules/nextflow/src/test/groovy/nextflow/script/ScriptMetaTest.groovy` (extend; create if absent)

**Interfaces:**
- Produces: `Set<String> ScriptMeta.getIncludedProcessNames()` — names of processes brought in via `include` (the `imports` map), excluding locally-defined processes.

### Steps
- [ ] **Write the failing test.** In `ScriptMetaTest.groovy`, follow the file's existing pattern for building a `ScriptMeta` with local definitions and imported components (look at how the existing `getProcessNames`/`getLocalProcessNames` tests construct fixtures — read the file first). Assert that for a meta with one local process `foo` and one imported process `bar`, `getIncludedProcessNames() == ['bar'] as Set` (only the import), while `getLocalProcessNames() == ['foo'] as Set`.
- [ ] **Run it — expect FAIL** (method missing):
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.script.ScriptMetaTest"
  ```
- [ ] **Add the method** in `ScriptMeta.groovy`, mirroring `getLocalProcessNames()` but iterating `imports` instead of `definitions`:
  ```groovy
      /**
       * @return the names of the processes brought into this script via
       * {@code include} statements (the imported components), excluding processes
       * defined locally in this script.
       */
      Set<String> getIncludedProcessNames() {
          final result = new HashSet<String>(imports.size())
          for( final item : imports.values() ) {
              if( item instanceof ProcessDef )
                  result.add(((ProcessDef) item).name)
          }
          return result
      }
  ```
  (If `imports`'s value type or `ProcessDef.name` access differs, match the exact style used by the adjacent `getProcessNames()` — read those lines first.)
- [ ] **Run it — expect PASS.**
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.script.ScriptMetaTest"
  ```
- [ ] **Commit.**
  ```bash
  git add modules/nextflow/src/main/groovy/nextflow/script/ScriptMeta.groovy modules/nextflow/src/test/groovy/nextflow/script/ScriptMetaTest.groovy
  git commit -s -m "feat: ScriptMeta.getIncludedProcessNames for include-driven agent tools"
  ```

---

## Task 2: `SandboxGuard` (path containment)

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/agent/SandboxGuard.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/agent/SandboxGuardTest.groovy`

**Interfaces:**
- Produces: `static boolean SandboxGuard.isAllowed(Path candidate, Path workDir, Collection<Path> readableDirs, boolean write)` — true iff `candidate`, after normalization/realpath, is inside `workDir` (for write) or inside `workDir` or any of `readableDirs` (for read); rejects `..` traversal, absolute paths outside the allowed roots, and symlink targets that escape.

### Steps
- [ ] **Write the failing tests** in `SandboxGuardTest.groovy`. Use JUnit/Spock `@TempDir` to create a real sandbox dir and a sibling "outside" dir. Cover: a file directly inside `workDir` (write=true → allowed; read=true → allowed); a path containing `..` that escapes (`workDir.resolve('../evil')` → rejected for both); an absolute path outside all roots (rejected); a path inside a `readableDirs` entry (read=true → allowed; write=true → rejected, since writes are workDir-only); a **symlink inside workDir whose target is outside** (create `Files.createSymbolicLink(workDir.resolve('link'), outsideFile)`; reading via the link → rejected). Full test:
  ```groovy
  /*
   * Copyright 2013-2026, Seqera Labs
   * Licensed under the Apache License, Version 2.0 (the "License"); ... (full Apache header)
   */
  package nextflow.agent

  import java.nio.file.Files
  import java.nio.file.Path
  import spock.lang.Specification
  import spock.lang.TempDir

  class SandboxGuardTest extends Specification {
      @TempDir Path tmp

      def 'should allow files inside the work dir for read and write'() {
          given:
          def work = Files.createDirectories(tmp.resolve('work'))
          def f = work.resolve('out.txt')
          expect:
          SandboxGuard.isAllowed(f, work, [] as Set, true)
          SandboxGuard.isAllowed(f, work, [] as Set, false)
      }

      def 'should reject parent-traversal escape'() {
          given:
          def work = Files.createDirectories(tmp.resolve('work'))
          def escape = work.resolve('../secret.txt')
          expect:
          !SandboxGuard.isAllowed(escape, work, [] as Set, false)
          !SandboxGuard.isAllowed(escape, work, [] as Set, true)
      }

      def 'should reject an absolute path outside all roots'() {
          given:
          def work = Files.createDirectories(tmp.resolve('work'))
          def outside = Files.createDirectories(tmp.resolve('outside')).resolve('x.txt')
          expect:
          !SandboxGuard.isAllowed(outside, work, [] as Set, false)
      }

      def 'should allow reads from a whitelisted module-output dir but not writes'() {
          given:
          def work = Files.createDirectories(tmp.resolve('work'))
          def mod = Files.createDirectories(tmp.resolve('moduleout'))
          def f = mod.resolve('result.fa')
          expect:
          SandboxGuard.isAllowed(f, work, [mod] as Set, false)
          !SandboxGuard.isAllowed(f, work, [mod] as Set, true)
      }

      def 'should reject a symlink whose target escapes the sandbox'() {
          given:
          def work = Files.createDirectories(tmp.resolve('work'))
          def outside = Files.createDirectories(tmp.resolve('outside'))
          def secret = Files.write(outside.resolve('secret.txt'), 'x'.bytes)
          def link = Files.createSymbolicLink(work.resolve('link.txt'), secret)
          expect:
          !SandboxGuard.isAllowed(link, work, [] as Set, false)
      }
  }
  ```
- [ ] **Run — expect FAIL** (class missing):
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.SandboxGuardTest"
  ```
- [ ] **Create `SandboxGuard.groovy`:**
  ```groovy
  /*
   * Copyright 2013-2026, Seqera Labs
   * Licensed under the Apache License, Version 2.0 (the "License"); ... (full Apache header)
   */
  package nextflow.agent

  import java.nio.file.Files
  import java.nio.file.Path

  import groovy.transform.CompileStatic

  /**
   * Pure path-containment checks for the agent {@code filesystem} tool. Writes are
   * confined to the agent work dir; reads are allowed within the work dir or any
   * of the per-invocation readable dirs (module-output dirs). Resolves symlinks and
   * normalizes so that {@code ..} traversal and symlink targets that escape the
   * sandbox are rejected.
   */
  @CompileStatic
  class SandboxGuard {

      static boolean isAllowed(Path candidate, Path workDir, Collection<Path> readableDirs, boolean write) {
          if( candidate == null || workDir == null )
              return false
          final real = realOf(candidate)
          final root = realOf(workDir)
          if( isInside(real, root) )
              return true
          if( write )
              return false
          for( final dir : (readableDirs ?: Collections.<Path>emptyList()) ) {
              if( isInside(real, realOf(dir)) )
                  return true
          }
          return false
      }

      /**
       * Real path of an existing target, or the normalized absolute path of the
       * nearest existing ancestor joined with the remaining (non-existent) tail —
       * so a not-yet-created write target is checked against its real parent (which
       * defeats symlink escape) rather than its literal lexical path.
       */
      private static Path realOf(Path p) {
          Path abs = p.toAbsolutePath().normalize()
          if( Files.exists(abs) )
              return abs.toRealPath()
          // walk up to the nearest existing ancestor, realpath it, re-append the tail
          Path existing = abs
          final tail = new ArrayList<String>()
          while( existing != null && !Files.exists(existing) ) {
              tail.add(0, existing.getFileName().toString())
              existing = existing.getParent()
          }
          if( existing == null )
              return abs
          Path base = existing.toRealPath()
          for( final seg : tail )
              base = base.resolve(seg)
          return base.normalize()
      }

      private static boolean isInside(Path child, Path root) {
          return child.startsWith(root)
      }
  }
  ```
- [ ] **Run — expect PASS.**
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.SandboxGuardTest"
  ```
- [ ] **Commit.**
  ```bash
  git add modules/nextflow/src/main/groovy/nextflow/agent/SandboxGuard.groovy modules/nextflow/src/test/groovy/nextflow/agent/SandboxGuardTest.groovy
  git commit -s -m "feat: SandboxGuard path containment for the agent filesystem tool"
  ```

---

## Task 3: per-invocation work-dir sandbox, `DispatchContext` + ThreadLocal, and the `filesystem` tool

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/agent/DispatchContext.groovy`
- Create: `modules/nextflow/src/main/groovy/nextflow/agent/FilesystemToolSchema.groovy`
- Modify: `modules/nextflow/src/main/groovy/nextflow/agent/ModuleToolBridge.groovy` (ThreadLocal context; `call()` routes `filesystem`; `callFilesystem`; capability-aware descriptor list)
- Modify: `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` (`createToolBridge` recognises the `'filesystem'` capability; `run()` allocates the work dir + sets/clears the ThreadLocal; cleanup)
- Test: `modules/nextflow/src/test/groovy/nextflow/agent/ModuleToolBridgeFilesystemTest.groovy`

**Interfaces:**
- `DispatchContext`: `Path workDir`; `Set<Path> readableDirs` (mutable; seeded with `workDir`).
- `FilesystemToolSchema.descriptor() : ToolDescriptor` — name `filesystem`, inputSchema `{path:string, operation:enum[read,write,list,exists], content:string?}`.
- `ModuleToolBridge`: `static void setContext(DispatchContext ctx)` / `static void clearContext()` (ThreadLocal); `call('filesystem', argsJson)` → `callFilesystem`.
- `AgentDef`: per-invocation work dir under `session.workDir`; sets the ThreadLocal around `runner.run`.

### Steps
- [ ] **Write the failing bridge test** `ModuleToolBridgeFilesystemTest.groovy`. Construct a `ModuleToolBridge` with NO modules but the `filesystem` capability enabled (see the constructor/factory you will add). Set a `DispatchContext` with a `@TempDir` work dir. Exercise `call('filesystem', json)` for: `write` (creates a file in the work dir → returns success JSON), `read` (returns the content), `list` (lists the dir), `exists`; and an out-of-sandbox `read` of `../x` → returns `{"error":...}`. Assert results. (Model the JSON arg shape on `FilesystemToolSchema` below.) Read an existing bridge test (e.g. `AgentToolBridgeIntegrationTest`) first to match construction/serialization idioms. Remember to `clearContext()` in `cleanup()`.
- [ ] **Run — expect FAIL.**
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.ModuleToolBridgeFilesystemTest"
  ```
- [ ] **Create `DispatchContext.groovy`:**
  ```groovy
  /* Apache header */
  package nextflow.agent

  import java.nio.file.Path
  import java.util.concurrent.ConcurrentHashMap
  import groovy.transform.CompileStatic

  /**
   * Per-agent-invocation dispatch context: the sandbox work dir and the set of
   * directories the filesystem tool may read (the work dir plus the output dirs of
   * modules run during this invocation). Created per input record by the agent
   * operator and threaded to {@link ModuleToolBridge} via a ThreadLocal, so the
   * shared, pre-ignition bridge holds no per-record state.
   */
  @CompileStatic
  class DispatchContext {
      final Path workDir
      final Set<Path> readableDirs

      DispatchContext(Path workDir) {
          this.workDir = workDir
          this.readableDirs = ConcurrentHashMap.newKeySet()
          if( workDir != null )
              this.readableDirs.add(workDir)
      }

      void addReadableDir(Path dir) {
          if( dir != null )
              readableDirs.add(dir)
      }
  }
  ```
- [ ] **Create `FilesystemToolSchema.groovy`:**
  ```groovy
  /* Apache header */
  package nextflow.agent

  import groovy.transform.CompileStatic

  /**
   * Descriptor for the generic {@code filesystem} agent tool. The work dir is bound
   * per-call from the dispatch context (never supplied by the LLM); the LLM provides
   * a path (relative to or inside the sandbox), an operation, and content for writes.
   */
  @CompileStatic
  class FilesystemToolSchema {

      static final String NAME = 'filesystem'

      static ToolDescriptor descriptor() {
          final input = [
              type: 'object',
              properties: [
                  path     : [type: 'string', description: 'File or directory path within the agent sandbox (the agent work dir) or a module-output path returned by a previous tool call.'],
                  operation: [type: 'string', enum: ['read','write','list','exists'], description: 'The filesystem operation to perform.'],
                  content  : [type: 'string', description: 'Content to write (required for the write operation; ignored otherwise).'],
              ],
              required: ['path','operation'],
              additionalProperties: false,
          ] as Map
          final desc = 'Read, write, list, or check files within the agent sandbox. Writes are confined to the agent work dir; reads may also target module-output paths returned by module_run. Small text files are returned inline; binary/large files are returned as path handles.'
          return new ToolDescriptor(NAME, desc, input, null)
      }
  }
  ```
- [ ] **Add the ThreadLocal context + `callFilesystem` + capability awareness to `ModuleToolBridge.groovy`.** Read the file first. Then:
  - Add a `private static final ThreadLocal<DispatchContext> CONTEXT = new ThreadLocal<>()` and `static void setContext(DispatchContext c){ CONTEXT.set(c) }` / `static void clearContext(){ CONTEXT.remove() }` / `private static DispatchContext context(){ CONTEXT.get() }`. Document that this relies on tool dispatch running on the calling thread (AiServices sequential; `executeToolsConcurrently` never enabled).
  - Add a `boolean filesystemEnabled` flag set via the constructor/factory; when true, `descriptors()` includes `FilesystemToolSchema.descriptor()`.
  - In `call(String toolName, String argsJson)` (the existing `synchronized` method), BEFORE the `tools.get(toolName)` lookup, branch: `if( toolName == FilesystemToolSchema.NAME ) return callFilesystem(parseArgs(toolName, argsJson))`. Keep the existing module-tool path for all other names. Keep the existing try/catch error wrapping.
  - Implement `private String callFilesystem(Map args)`: read `operation`/`path`/`content`; resolve the path against `context().workDir` if relative; gate via `SandboxGuard.isAllowed(resolved, ctx.workDir, ctx.readableDirs, write=(op=='write'))` → on failure return `JsonOutput.toJson([error: "path outside sandbox: ${path}"])`; if no context, return `{"error":"filesystem tool unavailable: no sandbox context"}`. For `read` use `ToolOutputReader.readOrHandle(resolved, maxInlineBytes)` and return `[content: ...]` (or the inline/path-handle shape); `write` writes `content` (UTF-8), returns `[path: <abs>, bytes: n]`; `list` returns `[entries: [names...]]`; `exists` returns `[exists: bool]`. If the work dir scheme is not `file`, return `{"error":...}`.
- [ ] **Wire the capability + work dir in `AgentDef.groovy`.** Read `createToolBridge` (lines 248-298) and `run()` (lines 135-214) first. Then:
  - In `createToolBridge`: split `declared` tools into capability strings vs legacy module refs. Recognise `'filesystem'` (→ set the bridge's `filesystemEnabled`). Build a bridge whenever ANY capability OR any legacy module is declared (so `tools 'filesystem'` alone returns a non-null bridge with no wired modules). Unknown capability-looking strings still fall through to the legacy module-resolution error path.
  - In `run()`: when the bridge has the `filesystem` capability, allocate a per-invocation work dir in the `mapper` closure: `final Path agentWorkDir = FileHelper.getWorkFolder(session.workDir, CacheHelper.hasher([session.uniqueId, name, idx, inputJson]).hash()); Files.createDirectories(agentWorkDir)`. Wrap the `runner.run(req)` call: `ModuleToolBridge.setContext(new DispatchContext(agentWorkDir)); try { result = runner.run(req) } finally { ModuleToolBridge.clearContext() }`. (`session` via `Global.session as Session`, as elsewhere in this file.)
  - In `afterStop`/bridge close: when `session.config.cleanup` is set, delete the per-invocation work dirs created for this agent. (Track created dirs on the AgentDef instance, or delete the agent-scoped bucket.) If cleanup tracking is non-trivial, leave the dirs (they live under `work/`) and note it — do not block the task on cleanup.
- [ ] **Run — expect PASS** (and compile the plugin/core):
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.ModuleToolBridgeFilesystemTest"
  ./gradlew :nextflow:compileGroovy
  ```
- [ ] **Boundary check:** `grep -RIl 'dev\.langchain4j' modules/nextflow/src/main` must be empty (the ThreadLocal/context are core-only). 
- [ ] **Commit.**
  ```bash
  git add modules/nextflow/src/main/groovy/nextflow/agent/DispatchContext.groovy \
          modules/nextflow/src/main/groovy/nextflow/agent/FilesystemToolSchema.groovy \
          modules/nextflow/src/main/groovy/nextflow/agent/ModuleToolBridge.groovy \
          modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy \
          modules/nextflow/src/test/groovy/nextflow/agent/ModuleToolBridgeFilesystemTest.groovy
  git commit -s -m "feat: agent work-dir sandbox + filesystem tool (ThreadLocal dispatch context)"
  ```

---

## Task 4: include-driven single `module_run` tool

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/agent/ModuleRunToolSchema.groovy`
- Modify: `modules/nextflow/src/main/groovy/nextflow/agent/ModuleToolBridge.groovy` (`module_run` descriptor; `callModuleRun` routes to the named internal module and records its output dirs)
- Modify: `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` (`createToolBridge`: `'module_run'` capability → enumerate `include`d modules, resolve + wire each, build the single `module_run` descriptor)
- Test: `modules/nextflow/src/test/groovy/nextflow/agent/AgentModuleRunToolTest.groovy`

**Interfaces:**
- `ModuleRunToolSchema.descriptor(List<String> moduleNames, Map<String,String> hints) : ToolDescriptor` — name `module_run`, inputSchema `{module: {type:string, enum:moduleNames}, args: {type:object, additionalProperties:true}}`, description aggregating per-module hints.
- `ModuleToolBridge.call('module_run', '{"module":"<name>","args":{...}}')` → looks up the internal `Tool` named `<name>`, runs it via the existing `callScalar`/`callSpec` with `args`, and adds each produced output-file parent dir to `context().readableDirs`.

### Steps
- [ ] **Write the failing test** `AgentModuleRunToolTest.groovy`. Use an in-scope typed process as the module (mirror `examples/agents/tool` / the existing in-scope bridge test). Build a bridge with `module_run` enabled over that one module. Assert: `descriptors()` contains exactly one tool named `module_run` whose inputSchema `module.enum` lists the module name; `call('module_run', '{"module":"uppercase","args":{"text":"hi"}}')` runs the process and returns its serialized output; `call('module_run', '{"module":"nope","args":{}}')` returns `{"error":...}` (unknown module). Read the existing in-scope bridge test first to reuse its process-wiring harness.
- [ ] **Run — expect FAIL.**
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.AgentModuleRunToolTest"
  ```
- [ ] **Create `ModuleRunToolSchema.groovy`:**
  ```groovy
  /* Apache header */
  package nextflow.agent

  import groovy.transform.CompileStatic

  /**
   * Builds the single {@code module_run} tool descriptor. The LLM picks a module by
   * name (constrained to an enum of the script's included modules) and supplies a
   * generic {@code args} object marshalled per that module's existing binding rules.
   * Per-module input hints are aggregated into the description so the LLM knows what
   * each module needs without a separate tool per module.
   */
  @CompileStatic
  class ModuleRunToolSchema {

      static final String NAME = 'module_run'

      static ToolDescriptor descriptor(List<String> moduleNames, Map<String,String> hints) {
          final input = [
              type: 'object',
              properties: [
                  module: [type: 'string', enum: new ArrayList<String>(moduleNames),
                           description: 'Which module to run. Must be one of the available modules.'],
                  args  : [type: 'object', additionalProperties: true,
                           description: 'The module inputs, as a flat object matching the chosen module (see the per-module input hints in this tool description).'],
              ],
              required: ['module','args'],
              additionalProperties: false,
          ] as Map
          final sb = new StringBuilder()
          sb.append('Run one of the available Nextflow modules as a real task and return its outputs as JSON. Available modules and their inputs:')
          for( final name : moduleNames ) {
              sb.append('\n\n### ').append(name).append('\n')
              sb.append(hints.get(name) ?: 'inputs: (see module documentation)')
          }
          sb.append('\n\nFile/path outputs are returned as absolute path strings (read them with the filesystem tool if needed), never file contents.')
          return new ToolDescriptor(NAME, sb.toString(), input, null)
      }
  }
  ```
- [ ] **Add `callModuleRun` + the `module_run` descriptor to `ModuleToolBridge.groovy`.** Read the file first. Then:
  - Track which internal `Tool`s are "module_run candidates" (all the wired modules) and a `boolean moduleRunEnabled` flag.
  - When `moduleRunEnabled`, `descriptors()` includes `ModuleRunToolSchema.descriptor(moduleNames, hints)` (build `hints` per module from `ModuleSpecToolSchema.inputSchema(spec)`+`outputDescription(spec)` or `ModuleMetadataToolSchema.description(metadata)` when available, else a compact rendering of `ProcessToolSchema.inputSchema(proc)` for a typed in-scope process). The individual per-module descriptors are NO LONGER added when `moduleRunEnabled`.
  - In `call(...)`, branch `if( toolName == ModuleRunToolSchema.NAME ) return callModuleRun(parseArgs(toolName, argsJson))`.
  - `private String callModuleRun(Map parsed)`: read `module` (String) and `args` (Map); `final tool = tools.get(module)`; if null → return `{"error":"unknown module: ${module} - available: ${moduleNames}"}`; else delegate to the existing `tool.spec != null ? callSpec(tool, args) : callScalar(tool, args)`. After the run, scan the produced output JSON for file path strings and add each file's parent dir to `context()?.addReadableDir(...)` so the filesystem tool may read module outputs. (Reuse the same path-detection ToolOutputReader uses, or record dirs inside `callSpec`/`callScalar` via the context — keep the bridge otherwise stateless.)
- [ ] **Wire the `module_run` capability in `AgentDef.createToolBridge`.** When `'module_run'` is among the declared capabilities: enumerate the included modules via `ScriptMeta.getIncludedProcessNames()` (Task 1), resolve each to its `ProcessDef` via `meta.getProcess(name)`, load any sibling `ModuleSpec` (reuse `loadSiblingSpec`), and wire them into the bridge exactly as the existing module path does (`wireSpec`/`wireScalar`) — but mark the bridge `moduleRunEnabled = true` so it advertises ONE `module_run` tool. (Legacy `tools 'nf-core/...'`/path/in-scope-name strings still resolve via the existing path and still advertise per-module tools — deprecated but unbroken.)
- [ ] **Run — expect PASS** + compile:
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.AgentModuleRunToolTest"
  ./gradlew :nextflow:compileGroovy
  ```
- [ ] **Commit.**
  ```bash
  git add modules/nextflow/src/main/groovy/nextflow/agent/ModuleRunToolSchema.groovy \
          modules/nextflow/src/main/groovy/nextflow/agent/ModuleToolBridge.groovy \
          modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy \
          modules/nextflow/src/test/groovy/nextflow/agent/AgentModuleRunToolTest.groovy
  git commit -s -m "feat: include-driven module_run tool (single tool, module enum)"
  ```

---

## Task 5: migrate examples + document the new model

**Files:**
- Modify: `examples/agents/skesa/main.nf`, `examples/agents/isolate-triage/main.nf`, `examples/agents/tool/main.nf` (and `two-agents` if it uses module tools)
- Create: `examples/agents/filesystem/` (a small example using `tools 'module_run', 'filesystem'`)
- Modify: `docs/agent.mdx` (Tools section, Directives table, Limitations)
- Test: `modules/nextflow/src/test/groovy/.../AgentScriptLoadingTest.groovy` — a script with `include {...}` + `tools 'module_run', 'filesystem'` loads

**Interfaces:** none (examples + docs).

### Steps
- [ ] **Migrate the example scripts** from the legacy `tools 'nf-core/skesa'` form to `include { skesa } from 'nf-core/skesa'` + `tools 'module_run'`. For `isolate-triage`, add `'filesystem'` and have the instruction mention writing a summary file. Keep the high-level instructions. Verify each migrated script PARSES/LOADS (a loading test or `nextflow inspect`/dry-run as available offline). Do not run real registry/container execution here.
- [ ] **Add the loading test** in `AgentScriptLoadingTest.groovy`: an agent declaring `include`d modules + `tools 'module_run','filesystem'` loads and the resolved `AgentDef.getTools()` contains the capability strings. Run:
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.script.parser.v2.AgentScriptLoadingTest"
  ```
- [ ] **Update `docs/agent.mdx`**: rewrite the "Tools (calling modules)" section for the capability model — `tools 'module_run', 'filesystem'`; modules come from `include` statements; the single `module_run(module, args)` tool; the `filesystem` tool (read/write/list/exists, sandboxed to the agent work dir + module-output paths); note `bash`/`nextflow run` are not available; note the legacy string form is deprecated. Update the Directives table `tools` row and the Limitations list.
- [ ] **Commit.**
  ```bash
  git add examples/agents docs/agent.mdx modules/nextflow/src/test/groovy/nextflow/script/parser/v2/AgentScriptLoadingTest.groovy
  git commit -s -m "docs+examples: migrate agents to capability tools (module_run + filesystem)"
  ```

---

## Task 6: Milestone 3 gate

**Files:** none modified — verification only.

### Exit-criteria steps
- [ ] **Boundary grep: core has no langchain4j.**
  ```bash
  grep -RIl 'dev\.langchain4j' modules/nextflow/src/main || echo "CLEAN: no langchain4j imports in core"
  ```
  Expected: `CLEAN: no langchain4j imports in core`
- [ ] **Boundary grep: plugin untouched (no ProcessDef/Channel/Path).**
  ```bash
  grep -RIn -E 'import nextflow\.(script\.ProcessDef|extension\.|Channel)|import java\.nio\.file\.Path' plugins/nf-agent/src/main || echo "CLEAN: plugin does not touch ProcessDef/Channel/Path"
  ```
  Expected: `CLEAN: ...`
- [ ] **Guard: executeToolsConcurrently never called.**
  ```bash
  grep -RIn 'executeToolsConcurrently' plugins/nf-agent/src/main || echo "CLEAN: executeToolsConcurrently never called"
  ```
  Expected: `CLEAN: ...`
- [ ] **New M3 unit tests.**
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.script.ScriptMetaTest" --tests "nextflow.agent.SandboxGuardTest" --tests "nextflow.agent.ModuleToolBridgeFilesystemTest" --tests "nextflow.agent.AgentModuleRunToolTest" --tests "nextflow.script.parser.v2.AgentScriptLoadingTest"
  ```
  Expected: `BUILD SUCCESSFUL`.
- [ ] **Full nextflow + plugin agent suites (regression — legacy module-string tests must still pass).**
  ```bash
  ./gradlew :nextflow:test --tests "nextflow.agent.*" && NXF_SMOKE=1 ./gradlew :plugins:nf-agent:test
  ```
  Expected: `BUILD SUCCESSFUL`, 0 failures.
- [ ] **Gated E2E (only if `OPENAI_API_KEY` is set).** Run an agent using `tools 'module_run', 'filesystem'` end-to-end if a gated test exists / can be added; otherwise confirm existing gated e2e still pass.
  ```bash
  if [ -n "$OPENAI_API_KEY" ]; then ./gradlew :plugins:nf-agent:test --tests "nextflow.agent.AgentToolEndToEndTest" --tests "nextflow.agent.AgentEndToEndTest"; else echo "SKIP: OPENAI_API_KEY not set"; fi
  ```
- [ ] **Clean worktree** (`.superpowers/` ignored; `examples/agents/two-agents/` pre-existing untracked unless migrated).
  ```bash
  git status --porcelain
  ```

### Milestone 3 exit criteria (all must hold)
1. `tools` accepts capability strings `'filesystem'`/`'module_run'`, validated at run time.
2. Module tools for `module_run` are sourced from `include` statements (`ScriptMeta.getIncludedProcessNames()`), surfaced as ONE `module_run(module, args)` tool with a `module` enum and aggregated per-module hints.
3. Each agent invocation gets a work-dir sandbox; the `filesystem` tool reads/writes within it (+ reads module-output dirs), enforced by `SandboxGuard`; per-record context is `ThreadLocal`, the bridge stateless across records.
4. Core imports no `dev.langchain4j.*`; the plugin is unchanged; `executeToolsConcurrently` never called.
5. `bash`/`nextflow run` deferred. Legacy module-string tools still work (deprecated); examples migrated; docs updated.
6. New unit tests green; full agent suites green; gated e2e green when credentials present.

## Open question for the user (non-blocking)
The legacy `tools '<module-ref-string>'` path is kept working-but-deprecated in M3 (removing it would churn ~7 core tests + the example tests). Confirm whether a follow-up should **fully remove** the legacy string-module path (and migrate those tests) or keep it as a permanent fallback.
