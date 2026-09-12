# Pluggable agent runners — the harness SPI and its implementations

- Status: implemented (experimental, unreleased)
- Scope: why the LLM harness is an extension point, the `AgentRunner` SPI, runner selection, the two
  shipped runners, and how the runner-native/brokered tool split falls out
- Companion documents: [agent design](agent-design.md) (the primitive itself),
  [agent RPC](agent-rpc.md) (the transport an out-of-process runner uses)

## 1. Why the harness is an extension point

`agent.runner` selects the LLM harness the way `process.executor` selects the executor. A pipeline
that hard-codes an HTTP call to one provider's endpoint has made that choice permanent; a pipeline
that declares an agent has not.

That matters beyond taste. The model endpoint being swappable is what lets the same pipeline target a
hosted API or a locally-hosted model without touching pipeline code — the local-first property the
agentic-genomics agenda asks for, delivered by the same mechanism that makes an executor swappable.

Two runners ship:

| Runner | Plugin | Loop runs | Container |
|---|---|---|---|
| `langchain4j` | `nf-agent` | driver JVM, on the dedicated `agent` executor | none |
| `pi` | `nf-agent-pi` | agent task container, on any executor | required |

## 2. The SPI

```groovy
interface AgentRunner {
    default String getName()                                   // stable: 'pi', 'langchain4j'
    default AgentLaunchSpec getLaunchSpec()                     // null => in-JVM runner
    default AgentRpcRegistration register(AgentRunnerRequest, boolean remote)
    String run(AgentRunnerRequest request)                      // in-JVM call path
}
```

Everything crossing this boundary is portable: `ToolDescriptor` and `SkillDescriptor` maps, an
`AgentRunnerRequest`, and a `String`/JSON result. Core never imports an LLM client; a runner never
touches `ProcessDef`, `Channel` or `Path`. `ToolDispatcher` in particular **never crosses** — the
host-side broker is the only component allowed to invoke it.

`getName()` is a defaulted method so existing SAM-shaped test runners stay source-compatible.
`getLaunchSpec()` returning `null` is the discriminator between the two execution shapes: a null spec
means the in-JVM path, a non-null one means the agent lowers to a `SCRIPTLET` submitted through the
configured executor.

### 2.1 Selection

`AgentRunnerProvider.get(name)` resolves an exact name. The no-argument form remains for tests and
single-runner compatibility but **fails on ambiguity instead of silently taking PF4J priority** —
priority is an unsafe user-facing selection mechanism once more than one runner is installed.

`AgentDef` resolves the runner once while building the task and folds the selected name into the
canonical task source, so **switching runners invalidates `-resume` entries**. Equal prompts sent
through different harnesses are not reproducibly equivalent, so runner identity belongs in the cache
key.

### 2.2 Where the request is assembled

Endpoint, credential and API provider are resolved **once, in core** (see
[agent design §7.1](agent-design.md#71-endpoint-and-credentials)) and handed to whichever runner is
selected. Credential resolution deliberately does **not** depend on runner selection: switching
`agent.runner` silently changing which environment variable is read would break the "resolved once in
core" invariant. Only when core resolves nothing may a runner resolve for itself — which `pi` can and
`langchain4j` cannot.

## 3. `langchain4j` (`nf-agent`)

The in-JVM runner, driven by langchain4j's **`AiServices`** API rather than a hand-rolled `ChatModel`
tool loop. `AiServices` owns the turn loop, tool advertisement, tool-result appending and the
structured-output response schema, so the plugin shrinks to model construction plus the adapters that
map portable descriptors onto langchain4j types:

- `ChatModelFactory` — `provider/model` → a langchain4j `ChatModel`, applying strict JSON schema when
  an output schema is present;
- `ModuleToolAdapter` — `ToolDescriptor` → `ToolSpecification`;
- `JsonSchemaMapper` — a portable schema `Map` → `JsonSchema`;
- `SkillAdapter` — `SkillDescriptor` → langchain4j `Skills`, providing `activate_skill` and
  `read_skill_resource` as an additive `ToolProvider`, with the skill catalog appended to the system
  message (the `skill_name` parameter is free text, not an enum — without the catalog the model cannot
  name a skill and activation never fires).

The runner supports the **OpenAI wire protocol only**. `openai/` names the protocol, so with
`agent.baseUrl` it reaches any endpoint that speaks it — vLLM, Ollama, llama.cpp, a gateway. Azure
OpenAI is *not* reachable this way: its URL shape and `api-version` handling need a client the plugin
does not ship.

**The runner must throw on failure, never return an error object.** Retries, `errorStrategy` and the
whole `resumeOrDie` path key on a `Throwable` reaching `task.error`; an error object would be bound as
a successful result. Provider-error retries should be capped low, since the OpenAI SDK already retries
internally and a 429 storm under fan-out compounds.

`shell:bash` is **rejected at agent-build time** on this runner, naming `pi`. The loop runs in the
driver JVM, so a bash tool would execute model-authored commands on the driver host with no container
boundary — and on a Kubernetes or Batch deployment the driver pod holds none of the pipeline's tooling
anyway. This is the one leaf where the cross-runner portability promise does not hold, and it is
deliberate.

## 4. `pi` (`nf-agent-pi`)

### 4.1 Why a child process

| Option | Verdict |
|---|---|
| Embed the JS runtime in the Nextflow JVM | Rejected — Pi targets Node APIs; Graal/JS embedding adds a second compatibility surface and weak dependency isolation |
| **Per-invocation child process with a bidirectional protocol** | **Selected** — natural Node runtime, live reverse tool calls, explicit cancellation and isolation |
| An ordinary file-in/file-out Nextflow process | Insufficient — a file-oriented task cannot synchronously call the live `ToolDispatcher`, and multi-turn state is lost |
| Turn-at-a-time canonical processes | Deferred — requires serializing and restoring full model/session state and scheduling a task per turn; high latency, provider-specific |
| Long-lived remote runner service plus callback gateway | Deferred — authentication, routing, tenancy isolation, persistence and failure recovery well beyond a runner SPI |

The selected protocol works identically for a local container and a remote one, and leaves a clean
path to a remote transport implementing the same messages.

### 4.2 The harness

`harness/runner.mjs` drives `@earendil-works/pi-coding-agent` through its public SDK:

- an in-memory session per invocation;
- a system-prompt override carrying `instruction`, optional `goal` and the skill catalog;
- **no built-in Pi coding tools** except the ones the agent's `tools` directive explicitly selected;
- one custom tool per brokered `ToolDescriptor`, which writes a `tool_call` frame and awaits its
  correlated `tool_result`;
- skill tools compatible with `activate_skill` / `read_skill_resource`;
- a terminating `final_answer` tool for structured output, carrying the Nextflow output schema.

**Structured output uses a terminating tool rather than a second structuring pass.** It avoids a lossy
re-read of the conversation and works uniformly with and without other tools. One corrective turn is
allowed when a model returns prose instead of calling it.

### 4.3 Packaging: the image *is* the distribution unit

The plugin originally vendored the whole runtime — 173 MB of `node_modules` (18,569 files) plus a
13 MB Go binary built for whichever architecture ran the build, unpacked at runtime with a
`deleteOnExit` hook per file. **The plugin artifact was ~186 MB**, against single-digit MB for a
typical Nextflow plugin, and releasing it required Go and npm toolchains.

Decision: **`nf-agent-pi` ships no runtime.** The proxy and harness live in a container image and the
jar is Groovy only. A container is therefore required for the `pi` runner on **every** executor,
including `local`.

That fixes all four defects at once — no vendored tree, no architecture matrix, no toolchain in the
release build, no extraction. Cross-compiling a `GOOS`/`GOARCH` matrix was rejected because it fixes
only the architecture problem and makes the size problem worse. Porting the Go proxy to Node is
attractive on its own terms and would delete a whole build stage, but it leaves the 173 MB, so it does
not make the plugin releasable; it remains worth doing later to reach one runtime.

**Containerization belongs to the executor layer, and that was achieved by deleting code.** The old
`createCanonicalBody` hand-rolled the decision between a driver-local and an in-container path. With
no host-local runtime there is one path, so the predicate, `shouldUseContainer`, `isOffloaded` and
`AgentLaunchSpec`'s local command pair all disappeared, and `BashWrapperBuilder` wraps the agent task
exactly as it wraps any other. Two bugs died with the predicate: `isOffloaded` treated *any* non-local
executor as containerized (wrong for a grid executor with no engine), and nothing consulted
`Executor.isContainerNative()`, so Kubernetes worked only incidentally.

What replaces the decision is **validation**. `AgentLaunchConditions.requireCanonicalLaunch` resolves
the `Executor` **instance** — not its name, because `isContainerNative()` can depend on the session
(the local executor is container-native under Fusion) and `containerConfigEngine()` decides which
engine block must be enabled — and rejects, before the run starts, a configuration that would not
containerize the task. The failure is a message naming `agent.container` and the engine, not a
`No such file` from inside the container.

| Agent executor | Required |
|---|---|
| `local`, or a grid executor | an enabled engine |
| container-native (`k8s`, `awsbatch`, …) | nothing further |

The image sets `ENTRYPOINT []` deliberately: canonical Nextflow/Fusion wrappers must retain entrypoint
control, so the generated agent script `exec`s the proxy explicitly.

### 4.4 The image is a release artifact, and the plugin declares it

If the image is the distribution unit, then leaving it outside the release makes the jar and its
runtime two deliverables that can drift — and they did: a proxy change landed after `VERSION` named a
tag that had already been published by hand, so the published `0.4.1` did not describe the tree.
Four decisions close that.

**The plugin declares the coordinate it needs.** A build-time generator reads the registry and image
name out of `build-image.sh` and the tag out of `plugins/nf-agent-pi/VERSION`, writing
`META-INF/nf-agent-pi-image.properties` into the jar; `PiAgentRunner` reads it back through a new
`AgentRunner.getDefaultContainer()` (defaulting to `null`, so `langchain4j` and every closure-coerced
test runner are untouched). What the jar asks for is by construction what the release publishes. A jar
without the resource yields `null` rather than throwing — degrading to core's existing
"must declare a container" error instead of an exception out of a static initializer inside PF4J
extension loading.

**`agent.container` becomes optional for `pi`.** `AgentDef` fills it from
`selectedRunner.getDefaultContainer()` when, and only when, the key is absent — after the
`executor`/`maxForks` defaults and before anything reads `config.get('container')`, so the
containerization guard, the per-task re-check and the container fingerprint all see one kind of value
and know nothing about the default. The test is `containsKey`, **not** truthiness:
`agent.container = false` keeps meaning "no container" and keeps failing.

> **Cache-key consequence.** The resolved container is not in `canonicalAgentSource` or
> `toolsFingerprint`; it reaches the hash generically, through `TaskHasher` adding
> `task.getContainerFingerprint()` for a containerized task. Existing agents are unaffected, because a
> canonical agent with no image is rejected *today* — every configuration that currently runs sets
> `agent.container` explicitly, so the new branch never fires for them. But bumping the plugin
> `VERSION` changes the default coordinate and therefore invalidates the cache of every agent relying
> on it. That is correct — a different runtime is a different task — and it is why the coordinate must
> never be a floating tag. Do **not** also add the container to the agent-specific fingerprints: it
> would be redundant and would move existing hashes.

**Publishing is idempotent, and verifies before it skips.** The tag is a runtime pin, hence immutable:
`push` publishes only when the tag is absent and exits 0 when it is present, mirroring the contract
`releasePluginToRegistryIfNotExists` already gives the plugin jars, so a partially failed release is
re-runnable. The multi-arch assertion runs on **both** paths — an earlier version checked it only
after a push, so a tag that existed but was single-arch was reported as "nothing to do" and the
release completed, git tag included, against an image that fails to pull on the other architecture. An
unreadable manifest ("absent", "unreachable" and "unauthenticated" are indistinguishable) falls
through to build-and-push, which is the authority.

**The release publishes it first.** `release.sh` runs the image build as step 1, ahead of everything
else, because it is the only step that reaches a third-party registry and because nothing before it is
undoable — the S3 upload, the plugin-registry entries and the git tag all are. Placed first, a failure
leaves nothing published.

**A drift guard makes "context changed ⇒ `VERSION` bumped" checkable.** It derives the guarded paths
from the `!` allowlist in `.dockerignore` rather than repeating them, since that file is already the
single statement of what enters the image — so a test-only proxy change does not fire, and a
hand-written `plugins/nf-agent-pi/**` would. Paths are prefixed with the plugin directory, because
`.dockerignore` is relative to the build context while a git pathspec is relative to the repo root;
unprefixed, the guard is a silent no-op. Not knowing the answer is always a warning, never a failure:
no git, no `VERSION`, and shallow clones all skip — the shallow case explicitly, because a graft
boundary is parentless, so a path-limited log names the *boundary* commit and the diff from there is a
false failure on a correct tree.

The registry is `public.cr.seqera.io/nextflow`, expressed as the single named constant everything
else derives from — the same registry and namespace the release already publishes `nextflow/nextflow`
to. That is deliberate: the release job's existing `public.cr.seqera.io` login authorizes the runner
image push too, so the image needs no credential of its own.

One operational note: the image build needs QEMU. The node stage is not `--platform=$BUILDPLATFORM`,
so `apt-get` and `npm ci` run on the target platform and the arm64 leg needs binfmt on an amd64
runner.

Delivering the image through Wave, building it on demand from the Dockerfile so each platform gets its
own architecture, remains the intended follow-up. The open question there is how Wave gets
build-context access to `agent-rpc/`, `package.json` and `package-lock.json`.

## 5. The tool split, per runner

**A tool reference is a contract, not an implementation.** Each runner satisfies it with whichever of
its own tools matches, and the wire name is the coordination point.

| Reference | Wire name | `langchain4j` | `pi` |
|---|---|---|---|
| `nf:module_run:X` | `X` | descriptor → `AiServices` tools; dispatch into `ModuleToolBridge` | brokered to the driver over RPC |
| `fs:read` / `write` / `ls` | `read` / `write` / `ls` | Groovy implementation in the driver JVM, behind `SandboxGuard` | SDK builtin, rooted at the session `cwd` |
| `fs:edit` / `grep` / `find` | `edit` / `grep` / `find` | Groovy implementation | SDK builtin |
| `shell:bash` | `bash` | **rejected at build time** | SDK builtin, inside the container |

The `pi` harness has **no local-vs-RPC branch**: every descriptor it receives is brokered,
unconditionally. Runner-native tools never become descriptors, so they cannot arrive that way — they
travel on their own request field, which **both** runners must consume (the `pi` harness passes the
names into the SDK allowlist; `langchain4j` rebuilds descriptors from them).

> **Known regression from that split.** The harness counts a tool turn inside each tool it *defines*,
> and the old mirrored filesystem tool was one of those. An SDK builtin is executed by the SDK and
> never passes through the counter, so on `pi` the whole `fs:` family and `shell:bash` are **not
> bounded by `maxIterations`**; the budget now covers only brokered module tools, skills and
> `final_answer`. The fix is to drive the counter off the session event stream
> (`tool_execution_start`), which counts every tool the model calls regardless of who executes it.

## 6. Credential delivery to a container

The resolved credential is delivered **in band** on the encrypted RPC start frame and installed as an
in-memory runtime key for that process only. It never enters the task environment, the task script or
the runner's credential store, and it takes precedence over the runner's own sources. Installing it
under the **model's** provider (the prefix the runner will dial), not under the driver-side resolved
`apiProvider`, is deliberate.

Sending is **gated on TLS**. With `agent.rpc.tls = false` the frame is cleartext, so a resolved
credential is withheld with a warning rather than shipped — `tls = false` stays a working escape hatch
and the out-of-band channels still reach the task.

This is what let the examples drop `agent.containerOptions = '-e OPENAI_API_KEY'` and their whole
`env { … }` blocks. Those channels remain documented for `tls = false` and for a provider whose
variable core does not resolve; the agent-scoped `containerOptions` form is preferable to `env`, since
a key in `env` also reaches every tool process.

## 7. Parity contract

Both runners are tested against the same behaviour contract: scalar and structured output; multiple
typed inputs and outputs; goal composition; module tools and recoverable dispatch errors; fatal tool
task failure; the filesystem sandbox; local skills and skill resources; tools plus structured output;
trace events and the resolved model; max-iteration failure; timeout and cancellation; parallel agent
execution without child-task starvation; and resume for tool-free/skills-only agents together with
resume opt-out for tool agents.

## 8. Residual risks

- **SDK API churn.** Pin exact versions and cover the harness with protocol tests.
- **Node/container startup overhead.** Measure before considering a pooled runner service — per-
  invocation isolation and correctness come first.
- **Provider differences in tool-schema support.** Fail capability checks before the first model call
  where the runner exposes enough metadata; otherwise surface a typed runner failure.
- **Secrets in child diagnostics.** Stderr is captured with a bounded size and redacted before
  inclusion in errors.
- **Cancellation while a module tool is running.** Session cancellation stays authoritative; the child
  is destroyed and the failure path clears interrupt state.
- **Which runner is the supported default** is still open. `pi` is the de-facto primary and every
  example selects it, while `langchain4j` has no launch spec and cannot offload at all. If
  `langchain4j` is the intended default *local* experience (in-JVM, no container), requiring a
  container engine locally for `pi` is uncontroversial; if `pi` is the default everywhere, that
  requirement deserves explicit sign-off.
