# Agent primitive — language surface and execution model

- Status: implemented (experimental, unreleased)
- Scope: the `agent` construct, its lowering to a task, typed I/O, tools, skills, modules,
  configuration, caching and lineage
- Companion documents: [agent runners](agent-runners.md) (the harness SPI and its
  implementations), [agent RPC](agent-rpc.md) (how an out-of-process runner reaches the driver)
- Decision record: [`adr/20260505-llm-agent-primitive.md`](../20260505-llm-agent-primitive.md)
- User guide: [`docs/agent.mdx`](../../docs/agent.mdx)

This document is the architectural reference for the primitive. It records *what was decided and
why*; the code is the reference for *how*. Where a decision has a non-obvious rationale or a known
sharp edge, that is stated here rather than left to be rediscovered.

## 1. The shape of the thing

An `agent` is a **task-shaped** primitive whose body is an LLM call, optionally with a tool-calling
loop. One input record in, one result out, executed as an ordinary `TaskRun`.

```nextflow
agent triage {
    model 'openai/gpt-5-mini'
    instruction 'You triage bacterial isolate assemblies step by step.'
    tools 'nf:module_run'

    input:  isolate: Isolate
    output: verdict: Verdict

    prompt: "Triage isolate ${isolate.sample_id} (${isolate.organism})."
}
```

Everything else in this document follows from one choice: **agent-as-task**. Because an agent
invocation is a `TaskRun` on the standard `TaskProcessor` machinery, parallelism, `maxForks`,
retries, `errorStrategy`, work directories, `-resume`, lineage, the progress table and the executor
abstraction are all *inherited* rather than rebuilt. The alternative considered and rejected — a
bespoke GPars operator with bolted-on progress and caching — would have had to hand-roll
`invokeTask → checkCachedOrLaunchTask → submitTask → monitor → finalizeTask → collectOutputsV2`,
strictly more code than reuse.

### 1.1 The boundary invariant

Two rules hold everywhere and are the reason the design stays composable:

- **Core never imports `dev.langchain4j.*`** (or any other LLM client). `modules/nextflow` and
  `modules/nf-lang` know only about the portable `AgentRunner` SPI: `ToolDescriptor` /
  `SkillDescriptor` maps, an `AgentRunnerRequest`, and a `String`/JSON result.
- **A runner plugin never touches `ProcessDef`, `Channel` or `Path`.** Dataflow, module resolution,
  schema derivation and JSON→record binding are core concerns.

Every mechanism below is designed to keep that seam intact. See [agent runners](agent-runners.md).

## 2. Lowering: agent → `TaskRun`

`AgentDef.run(Object[])` mirrors `ProcessDef.runV2` and produces the same three artifacts a v2
process produces:

1. a **`ProcessConfigV2`** with one `ProcessInput` per declared `AgentInput` and one `ProcessOutput`
   per declared `AgentOutput`;
2. a **`BodyDef`** — `exec` (Groovy, in-JVM runner) or `script` (canonical, out-of-process runner);
3. **invocation wiring** identical to `runV2`: one source channel per argument, `CH.create(singleton)`
   per output into a `LinkedHashMap`, then `createTaskProcessor().run()`, returning a `ChannelOut` so
   `myAgent.out.<name>` resolves.

`AgentDef` is deliberately **not** a `ProcessDef` subclass. It reproduces the ~12 lines of
`createTaskProcessor()`/`applyConfig()` instead, because making it a subclass would drag in
`ProcessDslV2`-driven config construction, which an agent cannot use: a process builds its
`ProcessConfigV2` once per definition, an agent builds one **per invocation** with synthesized
directives, body and output values.

### 2.1 Where the model call sits

| Runner kind | `BodyDef` type | Executed by | Model call runs |
|---|---|---|---|
| in-JVM (`langchain4j`) | `exec` → `ScriptType.GROOVY` | `AgentTaskHandler` on the dedicated `agent` executor | driver JVM |
| canonical (`pi`) | `script` → `SCRIPTLET` | whichever executor `agent.executor` names | agent task container |

The in-JVM body runs on a worker thread, **not** the GPars operator thread. That single fact is what
turns a serial operator into parallel tasks: the operator returns as soon as the task is submitted.

**The dedicated `agent` executor exists because of a two-throttle problem.** The stock
`LocalPollingMonitor` admits tasks only while `taskCpus <= availCpus`, so agent concurrency would be
capped at host CPU count — and would degrade to *serial* on a 1-vCPU CI runner, silently. Lifting
that requires **both** a non-CPU-gating monitor (`AgentPollingMonitor`) **and** a dedicated run pool
(`AgentTaskHandler`), because a non-gating monitor alone would admit handlers that then queue inside
the shared `session.getExecService()` pool — over-reporting "running" in the progress table while
only `poolSize` bodies actually execute. A dedicated pool also stops long LLM calls from starving
bash task launches, which hold a pool thread for their whole duration.

### 2.2 The output contract

`NativeTaskHandler`-style bodies put their return value in `task.stdout` — a single value, and one
an agent's declared output shares with its collected work-dir files. Neither can travel that way, so
**the body writes each output into `task.context`**, which `collectOutputsV2` resolves through
`TaskOutputResolver` and `bindOutputsV2` binds per channel.

The body closure runs `DELEGATE_ONLY` with `delegate = task.context`, and output names are dynamic
strings, so the write is `getDelegate().put(name, value)` — never bare assignment. `task.stdout` is
deliberately ignored on the in-JVM path. This is the highest-risk contract in the lowering and is
covered by focused tests.

## 3. Typed I/O

I/O uses the same typed declarations as a process and means the same thing. Exactly one deliberate
divergence exists (§3.3).

### 3.1 Inputs

N inputs, each `name: Type`. GPars' "fire when every input channel has a message" rule performs the
N-input combine natively, and `ProcessInputsDef.isSingleton()` decides fire-once vs fire-per-item —
so **fan-in needs no code at all**: `reducer(findings.collect())` feeds one value-channel item and
the task runs once.

`Path` inputs reach full process parity — staged into the task directory, bind-mounted in a
container, rendered as the plain filename in both the prompt and the input JSON. This works for a
bare `Path`, a `Set<Path>`/`List<Path>`, and `Path` fields of a record input, recursively.

**Path inference is shared at compile time, not reimplemented at runtime.** `_input_` receives a
`Class`, so `Set<Path>` arrives with its element type erased; runtime inference *cannot* reproduce
the `Collection<Path>` branch. `AgentToGroovyVisitor` therefore calls the same nf-lang helper
`ProcessToGroovyVisitorV2` uses, making agent staging identical to process staging by construction.

`tuple(...)` inputs are **rejected** with an explicit diagnostic rather than silently half-ignored,
and so is the destructured `record(...)` form.

### 3.2 Outputs

Two independent facts, and different consumers key on different ones:

- **An explicit right-hand side makes the output's value the RHS**, which removes it from the set the
  model answers. This is the process rule verbatim.
- **A `file(...)`/`files(...)` call in that RHS additionally registers an unstager**, making the value
  a work-directory collection.

| Declaration | Model answers it? | Collected from the work dir? |
|---|---|---|
| `answer: String` | yes | no |
| `report: Path = file('report.md')` | no | yes |
| `answer: String = 'fixed'` | no | no |

An agent whose outputs are *all* work-dir outputs is legal: the model is given no output contract and
its final text is discarded. **Nothing tells the model which filename to write**, so the prompt must
name it; a mismatch fails the task with an arity error. Automatic injection was rejected because it
would render an arbitrary pattern closure into instruction prose and change the resume cache key.

`AgentDsl` declares a nested `OutputDsl` exposing only `file`/`files`, so `eval()` and `stdout()`
fail to resolve in an agent output — a compile error, which is what we want. `ProcessUnstageVisitor`
gained a `filesOnly` mode for the same reason: its `env`/`eval` branches produce a `.command.env`
that only `BashWrapperBuilder` writes, which an in-JVM agent never has.

### 3.3 Structured output

An agent declares a single output — the model answers one value — so multiple values are returned
by declaring a `record` output. A named `record` output enables structured output: the record is
reflected into a JSON schema (`RecordSchema`) that constrains the model's response.

| Output | Schema |
|---|---|
| a scalar | none — free-text passthrough |
| a record | the **bare** `RecordSchema.of(type)` |

A record is **not** wrapped in an enclosing object. Wrapping it would burn one of OpenAI's five
nesting levels (dropping usable record depth 5→4) and change the wire shape, invalidating every
existing cache entry — for no benefit, since a record schema is already a valid strict object root.

Optional (`?`) record fields are handled by the runner, not by `RecordSchema`. Under OpenAI strict
mode langchain4j rewrites the field's type into the nullable union `["string","null"]` and forces
every field into `required` — the canonical strict-optional idiom. Encoding a union in
`RecordSchema` too would double-encode it and emit the more verbose `anyOf` shape, a *different*
wire schema and therefore a different cache key. This is pinned by a regression test rather than
fixed in source.

### 3.4 The prompt

`prompt:` is the templated user message. Like a process `script:` block it may contain helper
statements, and its **last expression** is what is sent. The whole input is additionally serialized
as JSON and appended, so the model sees structured values the template omits.

nf-lang captures the prompt closure's free-variable references into `PromptDef`, mirroring the
`VariableVisitor` that populates `BodyDef.valRefs` for a process body. Without that, a prompt closing
over a workflow global would be invisible to the cache key with no way to even warn.

## 4. Tools

### 4.1 The declaration grammar

`tools` entries are **namespaced references** — `family[:segment]*:name` — replacing the earlier flat
list that multiplexed capability literals, bare process names, module paths and registry coordinates
into one string and resolved them by trying each shape in turn.

| Family | Members | Executed by |
|---|---|---|
| `nf:` | `module_run`, itself a non-leaf over every in-scope process | the driver, as a real Nextflow task |
| `fs:` | `read`, `write`, `edit`, `ls`, `grep`, `find` | the runner |
| `shell:` | `bash` | the runner (`pi` only) |
| `mcp:` | reserved, unspecified | — |

Nine rules govern it; the load-bearing ones:

- **A non-leaf means its whole subtree**, so `nf:module_run` ≡ `nf:module_run:*`.
- **Globs only in the trailing segment, and only anchored to a family.** Bare `*` is rejected;
  `nf:*` and `fs:*` are legal because a family's membership is fixed by the Nextflow release rather
  than by remote configuration. A family with remote membership (`mcp:`) may never be globbed at its
  root.
- **Zero match is an error** — unknown family, missing explicit leaf, glob matching nothing, non-leaf
  over an empty subtree. A directive is a declaration, not a filter.
- **No negation operator.** The safety boundary sits on a family line instead: `shell:` escapes the
  filesystem sandbox, so it is its own family and `fs:*` never selects it. That is what makes "files
  but not a shell" expressible without `!`.
- **Order-independent set union**; overlapping references are idempotent.

The shape is Claude Code's permission-rule grammar with `:` for `__` and one deliberate deviation:
bare-prefix-implies-subtree is promoted from the server segment to *every* non-leaf depth. The colon
separator is not cosmetic — MCP tool names may contain `_`, so `__` is ambiguous, while `:` cannot
occur in a Nextflow identifier and is illegal in an OpenAI function name, forcing the
declaration/wire split below to be explicit.

### 4.2 Declaration namespace vs wire namespace

OpenAI function names are `[a-zA-Z0-9_-]{1,64}`, where a colon is illegal, so **a wire name never
contains a colon**. `nf:module_run:SAMTOOLS_SORT` is `SAMTOOLS_SORT` on the wire; `fs:read` is
`read`. The wire names are identical on both runners, so a pipeline stays portable.

Two checks run once at agent-build time, where skills and the output schema are also known:

- **Validate, never sanitize.** A selected process whose name is legal Nextflow but illegal on the
  wire (`my$proc`, over 64 characters) is a hard error. Silent renaming would collide `my$proc` with
  an existing `my_proc`.
- **One wire name from two different sources is a hard error naming both**, with a message stable
  under declaration order. The namespace being checked is the whole one the model sees, including the
  runner-injected `activate_skill`, `read_skill_resource` and `final_answer`.

### 4.3 Brokered vs runner-native

**`toolSpecs` carries brokered tools only.** This is the load-bearing invariant of the tool layer.

- **Brokered** (`nf:module_run:X`) — a `ToolDescriptor` crosses to the runner and the call comes back
  to the driver, which runs the module as a real Nextflow task wherever the agent runs.
- **Runner-native** (`fs:`, `shell:`) — resolved to *names* the runner is told to enable, carried on
  their own `AgentRunnerRequest` field, never as descriptors. The `pi` harness passes them to its SDK
  allowlist; `langchain4j` rebuilds descriptors from them and serves the tools in-JVM.

The alternative — a per-descriptor locality flag — was rejected: a mis-set flag would silently run a
container-side tool in the driver JVM, and it obliges the harness to reimplement what its runner
already ships.

**The sandbox boundary is per-runner and they are not the same boundary.** On `langchain4j`, `fs:`
tools go through `SandboxGuard`, confined to the task work dir for writes and additionally to
per-invocation whitelisted module-output paths for reads. On `pi` the SDK builtins are rooted at the
session `cwd` with the container as the outer bound — coarser — and `shell:bash` has no boundary
inside the container at all. This is documented rather than papered over.

### 4.4 The module tool bridge

A tool call must run synchronously from inside the agent body: the model is blocked waiting. But
invoking a `ProcessDef` after `session.fireDataflowNetwork()` **deadlocks** — a new operator's start
is deferred onto the igniter list, which is drained exactly once at ignition. Creating a dataflow
node post-ignition is not available to user code, and this is precisely why the primitive must live
in the engine rather than in a pipeline.

`ModuleToolBridge` is therefore built **in the workflow body, before ignition**. For each tool it
derives a `ToolDescriptor` and starts a dataflow **gateway** over a request queue. Each call carries
its own reply variable; the gateway creates fresh input/output channels and invokes a cloned
`ProcessDef`. Correlation is represented by dataflow variables rather than by ordering on a shared
lane, so **independent requests — including repeat calls to the same tool — execute concurrently**.
Request-scoped processors created after ignition start through the isolated
`Session.addProcessorIgniter` path, leaving workflow igniter behaviour unchanged.

> An earlier design pre-wired one persistent process instance per tool and made dispatch
> `synchronized`, because binding two calls' arguments before reading would break the 1:1
> input→output correlation on broadcast channels. The gateway replaced that. Serialization of tool
> calls today is a **runner** property — both harnesses execute tool calls sequentially — not a
> bridge property.

Descriptor sources, in order of preference:

| Source | When |
|---|---|
| registry `ModuleMetadata` | the module resolves in the registry — richest: per-field descriptions, patterns, enums, the nf-core `meta.id` convention |
| sibling `meta.yml` `ModuleSpec` | a local module, or registry metadata unreachable |
| declared typed process I/O (`ProcessToolSchema`) | an in-scope process with no `meta.yml` |

Schemas are **flattened**: `tuple(val(meta), path(fastq))` contributes one top-level property per
component, so the model passes `{"meta": …, "fastq": "/abs/path"}` rather than a nested array. A
drift guard warns when registry metadata and the installed `meta.yml` disagree on flattened input
names.

Per call: parse arguments; marshal them through `ProcessEntryHandler.getProcessArguments` (the same
path as `nextflow module run`, so type coercion and tuple assembly are shared); bind onto the
gateway; block on the reply; serialize outputs. `eval`/`topic`-routed bookkeeping outputs (nf-core
`versions`) are skipped authoritatively from the `ProcessDef`, because their channel never binds a
per-call value and reading it would block forever.

File outputs are **opaque absolute path strings** — the model chains them between tools and never
reads contents — with one deliberate exception: small structured text (`.json`/`.tsv`/`.csv`/…, under
`agent.maxToolOutputInlineSize`, non-binary) is **inlined**, which is what enables data-driven control
flow such as gating on a QC statistic.

Dispatch-level failures (unknown tool, unparseable arguments, bad marshalling) are returned as
`{"error": …}` so the model can correct and retry, counting toward `maxIterations`. A failing tool
**task** is not recovered — it surfaces through the process's own dataflow operator and aborts the
session. That asymmetry is a known limitation, not an oversight.

## 5. Skills

A skill is a `SKILL.md` folder — YAML frontmatter with `name` and `description`, then a Markdown body
— optionally carrying bundled resources.

- **Local**: a bare name resolves to `<definingDir>/skills/<name>/`. A single entry whose directory
  has no `SKILL.md` expands to every immediate subdirectory that has one.
- **Remote**: an explicit GitHub host (`github.com/<org>/<repo>[@rev]` and the `https://`/`git@`
  forms) is shallow-cloned into a **rev-keyed** cache directory. The bare `org/repo` form is
  deliberately *not* accepted — the registry-module rule already claims it.

Skills cross the SPI as portable `SkillDescriptor`/`SkillResource` DTOs with content eagerly loaded,
so no `Path` crosses and core stays LLM-client-free. The model sees each skill's name and description
up front, reads the body through `activate_skill`, and resources through `read_skill_resource`.
**Skills run no code**; the clone executes nothing from the repository.

Resource scanning skips `.git/` and symlinks, rejects paths escaping the skill directory, and caps at
64 files / 256 KB. Frontmatter is split by hand rather than through `fromYaml`, tolerating a BOM,
blank lines and CRLF so a skill that loads in langchain4j's own loader also loads here.

> A remote `SKILL.md` becomes model instructions — a prompt-injection surface, worsened by a moving
> ref. Pin a commit SHA.

## 6. Agents as modules

An agent is authorable as a module and consumed with the ordinary `include` statement. The headline
finding of that work: **the include path already worked**, because it is generic over `ComponentDef`
and never casts to `ProcessDef`, and the language layer already enumerated agents as includable
definitions. The deliverable was tests, docs and a handful of small parity fixes, each of which
*deleted an asymmetry* between an agent and a process — a stable `baseName` so selectors match an
aliased agent, agent names reaching the selector registry, an nf-lang formatter that was **deleting**
agent declarations, and a compile-time arity check.

The design principle: **a module's consumer cannot edit the module file.** Everything the consumer
needs from outside — aliasing, config selectors, comprehensible errors — must work without touching
it. That is the only reason any code changed.

Locked semantics:

- **File-relative resolution uses the *defining* module directory** — skills, relative tool paths,
  `moduleDir` — including under an alias, because `cloneWithName` preserves `owner`. There is **no
  fallback** to the project directory: a fallback would make a module's behaviour depend on the
  including project's layout, and would let a consumer silently shadow a module's instructions.
- **`nf:module_run` is module-lexical.** An agent sees only processes defined or included by the
  script that declares it. If the tool surface depended on the includer, the tool schema — and
  therefore the model's behaviour and the cache content — would vary by consumer, which is not a
  module. It is also the same rule a process body already obeys. Consequence to state loudly: **a
  module agent must include its own tool modules.**
- **Params are inherited** from the run; `params()`/`addParams()` are dead on the v2 path.
- Registry-hosted agent modules are designed but deferred: the transport works, publishing and
  validation do not.

`AgentDef.cloneWithName` is a bare `Object.clone()`, so aliases **share** directives, inputs, outputs
and the prompt with the template. This is safe because those fields are read-only after construction
and `buildAgentTask` runs single-threaded during graph construction. Do not "fix" it into a deep copy
without a failing test.

## 7. Configuration

The `agent` scope has two axes — **agent options** and **task directives** — and is fully independent
of `process`: an agent never reads `process` defaults or selectors, and vice versa.

Parity with the process scope (selectors, regex and `!` negation, `withLabel:`, `ext` merge,
repeatable and dynamic directives, selectors inside `profiles`) is achieved by **reusing
`ProcessConfigBuilder` verbatim** plus the smallest possible glue: a scope-name list, a reflective
agent-only key set, two inert constructor parameters and a `kind` noun for user-visible messages. No
second selector engine, no copy of `matchesSelector`/`ext`-merge/`putWithRepeat`.

The agent-only key set is **derived reflectively from `AgentConfig`**, not hand-listed, so an option
added to the scope cannot drift into being treated as a bogus task directive.

Two consequences worth stating:

- The `agent` scope uses the option names `AgentConfig` declares (`model`, `maxIterations`), not the
  body-directive names. Writing `agent { instruction = … }` is a user error and is reported as one,
  the same diagnostic a process typo gets. No alias layer was introduced.
- `agent.executor` defaults to `local`, **never** to the global `executor.name`. An in-JVM runner
  accepts no other value and runs on the dedicated `agent` executor (§2.1). `agent.container` falls
  back to the runner's own declared image when the key is absent — see
  [agent runners §4.4](agent-runners.md#44-the-image-is-a-release-artifact-and-the-plugin-declares-it).
- `agent.rpc.*` is read once per session, so a value inside a selector block resolves but has no
  effect.

### 7.1 Endpoint and credentials

One defect drove this design: `prefix == 'openai'` gated credential resolution, endpoint resolution
*and* client selection. One prefix answered three questions — which **wire protocol** to speak, whose
**credential** this is, whose **endpoint** this is. Fusing (1) with (2)–(3) is why Anthropic and
OpenRouter had no expressible spelling.

`agent.apiProvider` splits them out: it names the **credential and endpoint namespace** and does not
touch protocol selection, which the model-id prefix keeps.

```
provider = agent.apiProvider              // explicit always wins
        ?: inferFrom(neutral baseUrl)     // exact host or dot-suffix, never substring
        ?: prefixOf(model)
```

Inference outranks the prefix because the prefix names a *protocol*, so it is the weaker signal for
"whose key". Circularity is broken by inferring only from the **provider-neutral** base URL.

Both `apiKey` and `baseUrl` then resolve on the same four tiers: config → `NXF_AGENT_*` →
`<PROVIDER>_*` → unset. Tiers 1 and 2 are provider-neutral and reach any endpoint. **Tier 3 is
gated**, because a runner installs the credential it is handed ahead of anything it could resolve
itself — so an unrestricted tier 3 would be a credential-disclosure primitive. A `<PROVIDER>_API_KEY`
travels only when the endpoint agrees: it is that provider's host, it came from the same namespace's
`<PROVIDER>_BASE_URL`, `agent.apiProvider` is explicit and the host is not a *different* known
provider's, or no endpoint resolved at all **and** the provider equals the model prefix's own.

"Withheld" is not "missing". Both used to answer `null` and both used to become
`PLACEHOLDER_API_KEY`, converting a diagnosable misconfiguration into an opaque 401.
`credentialWithheld` carries the distinction, and the placeholder now exists only for a genuine
no-credential local endpoint. The two runners then diverge deliberately: `langchain4j` **errors**
(core is its only credential source), `pi` **warns once** (it has a store and a container environment
core cannot see, and a driver-side abort would break a working deployment).

## 8. Caching, resume and lineage

### 8.1 The cache key

No `TaskHasher`, `HashBuilder` or `CacheDB` edits. `TaskHasher.compute()` already hashes
`task.source` and the resolved per-record `task.inputs`, so **the single injection seam for static
agent identity is `BodyDef.source`**, set to a canonical string. It folds:

runner · resolved model id · resolved endpoint · instruction · goal · `maxIterations` · prompt
template and its captured variables · output schema · a name-sorted content fingerprint of every skill
· a fingerprint of every module tool (schema + backing process script) · the resolved runner-native
tool refs plus the runner's identity and version.

Three notes on that list:

- It folds **no agent name and no file path**, so moving or renaming a module directory does not
  invalidate; aliasing does, via the fully-qualified `task.processor.name` in `TaskHasher`.
- **Runner-native tools need their own contribution** precisely because they have no descriptor.
  Without it, an agent's `fs:*` would be absent from its key and a runner upgrade that changed what
  `edit` does would replay a stale generation. Runner identity + version is deliberately coarser than
  a schema hash — it is the granularity at which those tools' behaviour actually changes — and it
  correctly makes the key runner-dependent, since the same `fs:read` is a different implementation on
  each runner.
- **Globs are forward-open.** `nf:*` silently acquires tools added in a later release. The grammar
  does not mitigate that, but the key covers the **expanded** set, so gaining a tool re-runs the agent
  rather than replaying. Explicit refs are the reproducible choice.

Resume required exactly one shared-code change: `TaskRun.hasCacheableValues()` had to persist an
`exec` task's context. Everything else is inherited. `cache false` is the documented opt-out and is
free (`isCacheable()` false → no lookup, no storage).

**Resume replays a stored generation; it does not make the generation reproducible.** A floating model
alias warns on a cache-writing run, but pinning a dated snapshot is advisory, not enforced. Whether
the engine should ever refuse an unpinned model is an open question, and it is the closest thing here
to a validation-tier boundary.

Temperature is pinned to 0 for reproducible replays. For a `collect()`-fed reduce, `ArrayBag`
implements `Bag` and `HashBuilder` hashes it order-independently, so the key is stable across fan-in
reordering — but the *prompt* item order can still vary, so `toSortedList()` is the recommendation
when reproducible *fresh* runs matter.

### 8.2 Lineage

Zero new observers. A genuine `TaskProcessor`/`TaskHandler`/`TaskRun` routed through the monitor
fires the same `Session.notifyTask*` events that `WorkflowStatsObserver` and `LinObserver` already
subscribe to. An agent emits an **`AgentRun`** record naming the runner, the requested and resolved
model, instruction, goal, prompt template, iteration ceiling, output schema, and the resolved **wire**
tool names.

Recording the *declared* refs alongside the resolved ones was considered and rejected: it needs a new
field threaded through three layers for no consumer that exists today, and the expanded set is the one
that determines behaviour. Consequence: `fs:*` and `fs:read, fs:write, …` leave identical traces.

Not recorded: the rendered prompt (the template plus the recorded inputs reconstruct it), the resolved
command (it embeds the RPC capability token), token/turn/tool-call counts (not instrumented), and
`resolvedModel` for an out-of-process agent.

## 9. Known limitations and open surfaces

- **Tool-task failure aborts the run.** Only dispatch-level errors are recoverable by the model.
- **Tool calls are serialized per agent** by both runners; parallel/correlated tool calls within one
  invocation are a correlation question, deferred.
- **`nf:module_run` exposes declared process inputs only.** A stock nf-core module hides its tuning
  flags in `task.ext.args`, which is not in the schema, so an agent cannot tune it without a wrapper
  process that surfaces the knob as an input.
- **In-JVM agents do not materialize staged files.** Staging is `SCRIPTLET`-only; an `exec` body goes
  to a handler with no wrapper and no stage-in.
- **Agents whose inputs cannot be captured are pinned non-cacheable**, which is honestly
  non-replayable rather than falsely replayable. The boundary between keyable and non-keyable tools is
  a design surface, not a settled rule.
- **Validation is enabled, not obliged.** An agent result is an ordinary typed channel value, so it can
  flow into assertions, `nf-test` and QC processes. Nothing requires that it does.
- No duplicate-invocation guard; an agent-only module is not directly runnable; `nextflow inspect` does
  not see agents.

Registry- and community-side items — validation-tier declaration and enforcement, benchmarking
infrastructure, OOD/bias detector libraries, population and equity metadata, BCO bundle emission,
preregistration diffing — are metadata and governance work that this language change is meant to
enable, not absorb. See the ADR.
