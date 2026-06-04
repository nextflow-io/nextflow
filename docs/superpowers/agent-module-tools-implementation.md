# Agent Module Tools — Implementation Summary

How the Nextflow **agent** primitive exposes modules/processes as LLM **tools** that
execute as real dataflow nodes: how the agent's input reaches a specific module, how the
LLM drives tool invocation, and how a tool's outputs are wired back to the LLM and chained
into downstream tool executions.

> Scope: the **module-as-tool** execution path (spec-driven tools). Branch:
> `agent-syntax-and-model`. File/line references are approximate (use the method names).

---

## 1. The shape of an agent

```groovy
agent triage {
    model 'openai/gpt-5-mini'
    instruction '''high-level goal; the LLM decides which tools to call and when'''
    tools 'nf-core/skesa', 'nf-core/assemblyscan', 'nf-core/prokka'   // module tools
    input:  isolate: Isolate          // a typed record (the agent's input)
    output: verdict: String           // a plain value when tools are used
    prompt: """Triage '${isolate.sample_id}'. Reads: ${isolate.reads}"""
}
```

An `agent` is a dataflow operator like a process: it consumes an input channel and emits an
output channel. For each input item it runs **one LLM conversation**; inside that
conversation the LLM may call the declared `tools` any number of times. Each tool call runs
the corresponding **module as a real Nextflow task** (its container, work dir, caching) and
returns a result the LLM reads.

### Key files

| Concern | File |
|---|---|
| Agent DSL node, operator, tool resolution | `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` |
| Tool execution bridge (pre-wire + dispatch + output serialization) | `modules/nextflow/src/main/groovy/nextflow/agent/ModuleToolBridge.groovy` |
| Input marshalling (shared with `module run`) | `modules/nextflow/src/main/groovy/nextflow/script/ProcessEntryHandler.groovy` |
| Output content inlining | `modules/nextflow/src/main/groovy/nextflow/agent/ToolOutputReader.groovy` |
| Tool descriptor/schema from registry metadata | `ModuleMetadataToolSchema.groovy`, `plugins/nf-agent/.../JsonSchemaMapper.groovy` |
| LLM tool-calling loop | `plugins/nf-agent/src/main/nextflow/agent/LangChainAgentRunner.groovy` |
| Agent config scope | `modules/nextflow/src/main/groovy/nextflow/agent/AgentConfig.groovy` |

---

## 2. End-to-end data flow

```
input record ──► AgentDef MapOp ──► prompt + inputJson ──► AgentRunner (LLM loop)
                                                              │
                          ┌───────────────────────────────────┘
                          ▼ (per tool call)
        LLM emits tool name + flattened JSON args
                          │  ModuleToolBridge.call(name, argsJson)   [synchronized]
                          ▼
        ProcessEntryHandler.getProcessArguments  ──► [ [meta], file(reads) ]  (channel value)
                          │  bind onto pre-wired input queue(s)
                          ▼
        REAL module task runs on the executor (container/work dir/cache)
                          │  read pre-captured output channels (.val), skip eval/topic
                          ▼
        serialize outputs ──► JSON result   (small structured file => INLINED contents;
                          │                   bulk/binary file => absolute-path HANDLE)
                          ▼
        result fed back to the LLM as a ToolExecutionResultMessage
                          │  (LLM may chain a handle into the next tool, or reason over inlined data)
                          ▼
        LLM returns a final text answer ──► AgentDef emits it on the output channel
                          │
        agent input exhausted ──► bridge.close() poisons input queues ──► operators terminate
```

---

## 3. Wiring the tools (build time, **before** the dataflow network ignites)

Tools are resolved and pre-wired in the workflow body, before ignition. This is essential:
invoking a process synchronously *after* ignition deadlocks GPars; pre-wiring + a blocking
pull on the outputs is what makes synchronous tool calls work.

**`AgentDef.run()`** (`AgentDef.groovy`, ~L150–200):
- builds the prompt **mapper** closure and reads the `agent` config scope (model,
  `maxIterations`, `requestTimeout`, `maxToolOutputInlineSize`);
- calls **`createToolBridge()`** (~L220) which resolves every `tools` entry:
  - an **in-scope process name** → used directly;
  - a **local `.nf` path** or a **registry ref** (`nf-core/skesa`) →
    **`resolveRegistryTool()`** (~L412) / **`compileModuleProcess()`** (~L346) compile the
    module's `main.nf` into a `ProcessDef`. The public registry `ModuleMetadata`
    (`GET /api/v1/modules/{name}`) is fetched for the tool descriptor;
- wires the operator with **`applyAgentOperator()`** (~L511) and, on the operator's
  `afterStop`, calls **`bridge.close()`** to poison the tool input queues.

> **Critical (per-tool module isolation):** `compileModuleProcess` gives each tool module its
> **own** `ScriptBinding` (seeded with the owner's params) via `ScriptLoaderV2.setMainBinding`,
> mirroring `IncludeDef.loadModuleV2`. Without this, all tool modules shared `session.binding`
> and `BaseScript.setup()`'s `moduleDir` write was last-wins, so every tool's
> `conda "${moduleDir}/environment.yml"` resolved to the last-compiled module → one shared
> container for all tools.

**`ModuleToolBridge.wireSpec()`** (`ModuleToolBridge.groovy`, ~L209–270) — the heart of the
pre-wire, **once per tool**:
```groovy
// one persistent queue per declared input CHANNEL (== ProcessEntryHandler's arg count)
final toolIns = (0..<nInputs).collect { CH.queue() }
final cloned  = proc.cloneWithName(name)
final ChannelOut out = cloned.run(toolIns as Object[])      // pre-wire into the network

// CRITICAL: capture the output READ channels NOW, at wiring time
final outReadChannels = (0..<out.size()).collect { CH.getReadChannel(out[it]) }
```
The **descriptor** (tool name + description + input JSON schema) is built from the registry
`ModuleMetadata` via `ModuleMetadataToolSchema` (falling back to the sibling `meta.yml`
`ModuleSpec` offline). `JsonSchemaMapper` carries per-field descriptions, patterns and enums
into the langchain4j schema, and emits the nf-core `meta`→`{id}` sub-schema so the LLM can
build `meta.id` from the schema alone.

> **Critical (read-channel timing):** the output read channels are captured at **wiring**,
> not lazily at dispatch. These process outputs are broadcast-style — an adapter created
> *after* the task binds its value never receives it. Subscribing at wiring (before any value
> binds) is what lets `callSpec` read every output regardless of order; the scalar path
> (`wireScalar`, ~L203) always did this. Lazily subscribing per-read deadlocked the 2nd+
> output of multi-output tools (e.g. prokka).

---

## 4. Agent input → specific module input (marshalling)

### 4a. The agent's input becomes the prompt
`AgentDef.run()`'s **mapper** closure (~L182–194) runs once per input item:
```groovy
cl.setDelegate([(inputName): item]); final promptText = cl.call()?.toString()
final inputJson = toJson(item)                       // the input record as JSON
final req = new AgentRunnerRequest(model, instruction, promptText, maxIter,
                                   tools, outputSchema, inputJson, toolSpecs, bridge, timeout)
final result = runner.run(req)
```
The input record is rendered into the `prompt` **and** appended verbatim as JSON
(`composeMessages`, `LangChainAgentRunner.groovy` ~L145). The agent's input shape and a
tool's input shape are **independent** — the LLM bridges them.

### 4b. The LLM synthesizes the tool arguments
The LLM emits a tool call with a **flattened JSON** argument object shaped by the tool's
input schema, e.g.:
```
nf_core_skesa({"meta":{"id":"isolate_001"}, "fastq":"…/sample.fastq"})
```
It maps `isolate.sample_id → meta.id` and `isolate.reads → fastq` itself — no field-by-field
correspondence is required.

### 4c. The bridge marshals the flattened args onto the module's input channels
**`ModuleToolBridge.callSpec()`** (`ModuleToolBridge.groovy`, ~L431–445) reuses the **exact
same binding logic as `nextflow module run`**:
```groovy
final List args = ProcessEntryHandler.getProcessArguments(tool.processDef, parsed, spec)
for( int i=0; i<tool.toolIns.size(); i++ )
    tool.toolIns[i].bind(args[i])           // bind onto the pre-wired queues → task runs
```
**`ProcessEntryHandler.getProcessArguments`** (`ProcessEntryHandler.groovy`, ~L184/199 →
`getProcessArgumentsV1` ~L251): the flattened arg map is treated as the params map (dot-notation
folded into nested maps); each declared input is looked up by name, **type-coerced** from the
sibling `meta.yml` (`file`→`Nextflow.file`, `map`→Map, `integer`→Integer, …), and a tuple input
is assembled in order into a channel value such as `[[id:'s1'], file(reads)]`.

> **Critical (optional path inputs):** `getValueForInputV1` (~L341) treats a **null OR
> blank-string** file argument as *not provided* → empty list `[]` (the nf-core convention for
> optional `path` inputs). This avoids `file("")` throwing when the LLM passes `""` for an
> optional input it doesn't have (e.g. prokka's `proteins`/`prodigal_tf`).

A dispatch-level failure (unknown tool, unparseable JSON, bad argument) is **returned as a
`{"error": "..."}` tool result** rather than thrown (`ModuleToolBridge.call`, ~L388), so the LLM
can correct itself and retry instead of aborting the run.

---

## 5. LLM tool invocation loop

**`LangChainAgentRunner.runWithTools()`** (`LangChainAgentRunner.groovy`, ~L104–138) — the manual
tool-call loop (langchain4j; tools and structured-output are mutually exclusive in v1):
```groovy
final specs = request.toolSpecs.collect { ModuleToolAdapter.toToolSpecification(it) }
for( int i = 0; i < maxIterations; i++ ) {
    final response = model.chat(ChatRequest.builder().messages(messages)
                                   .toolSpecifications(specs).build())
    final ai = response.aiMessage()
    if( !ai.hasToolExecutionRequests() )
        return ai.text()                                   // final answer
    messages.add(ai)                                       // assistant turn (the tool requests)
    for( ToolExecutionRequest ter : ai.toolExecutionRequests() ) {
        final String result = request.dispatch.call(ter.name(), ter.arguments())   // ← the bridge
        messages.add(ToolExecutionResultMessage.from(ter, result))                 // feed result back
    }
}
```
- the tool specs are advertised on **every** chat request;
- `request.dispatch` is the `ModuleToolBridge` (a `ToolDispatcher`); `call(...)` is
  **`synchronized`** so concurrent tool requests in one turn stay input/output correlated;
- the loop ends when the model returns text with no tool requests, or hits `maxIterations`.

---

## 6. Tool output → result, and chaining into downstream tools

After the module task completes, **`callSpec`** reads the outputs and serializes them
(`ModuleToolBridge.groovy`, ~L446–474):
```groovy
final outParams = declaredOutputParams(tool.processDef)   // authoritative param list
for( int i=0; i<outputs.size(); i++ ) {
    final param = outputs[i]
    if( isEvalOutput(param) || isTopicOutput(outParams, i) ) continue   // skip versions/topic
    if( i >= nOut ) continue
    final ch = tool.outReadChannels[i]                    // read channel captured at WIRING
    final value = ch.val                                  // blocking pull
    result.put(key, serializeOutput(param, value))
}
return JsonOutput.toJson(result)
```

> **Critical (skip non-binding bookkeeping outputs):** an nf-core `versions` output is routed to
> a `topic:` (a topic-source channel that never binds a readable per-task value). It is detected
> **authoritatively from the ProcessDef** — `isTopicOutput` checks `OutParam.getChannelTopicName()`
> (guarded to `BaseOutParam`, since typed/V2 `ProcessOutput.getChannelTopicName()` throws) — not
> from the meta.yml `type`, which the registry types inconsistently (`eval` vs `string`). Reading
> such a channel's `.val` would block the dispatch forever.

### The two kinds of output — handle vs inlined
`serializeOutput` turns a tuple into a record keyed by component name; each file value goes
through **`ToolOutputReader.readOrHandle(path, maxInlineBytes)`** (`ToolOutputReader.groovy`):

- **bulk / binary data** (`.fa`, `.bam`, `.gz`, … — anything not in the text set, or over the
  size cap, or with NUL bytes) → returned as an **absolute path string** (an *opaque handle*).
  The LLM never sees the bytes; it forwards the handle to the next tool.
- **small structured text** (`.json`/`.tsv`/`.csv`/`.txt`/`.yaml`/… under
  `agent.maxToolOutputInlineSize`, default 32 KB) → its **contents are inlined** so the LLM can
  reason over them.

This split is what enables **downstream chaining** and **data-driven control flow**, e.g. in
`isolate-triage`:
1. `nf_core_skesa` → `fasta` is a `.fa` → returned as a **path handle**.
2. The LLM passes that handle as `nf_core_assemblyscan`'s `assembly` input (chaining).
3. `nf_core_assemblyscan` → its stats `.tsv` is small/structured → **inlined**; the LLM reads
   `n50_contig_length` / `total_contig` and applies the QC gate.
4. On PASS it calls `nf_core_prokka` (forwarding the contigs handle) and reports the `.gff`
   annotation path; on FAIL it stops and reports — prokka is never invoked.

The serialized JSON result becomes the `ToolExecutionResultMessage` content the LLM reads
next iteration.

---

## 7. Completion & termination

The agent operator runs one LLM conversation per input item. When the agent's input source is
exhausted, the operator's `afterStop` listener calls **`ModuleToolBridge.close()`**, which binds a
`PoisonPill` to every pre-wired tool input queue. The tool process operators receive the poison,
terminate, and propagate poison to their outputs (`TaskProcessor.sendPoisonPill`), so the dataflow
network unwinds and the session completes cleanly.

---

## 8. Critical code sections (quick index)

| What | Where | Why it matters |
|---|---|---|
| Pre-wire tools + **capture output read channels at wiring** | `ModuleToolBridge.wireSpec` | late per-read subscription deadlocks multi-output tools |
| Marshal LLM args → channel values (shared with `module run`) | `ProcessEntryHandler.getProcessArguments` / `getValueForInputV1` | flattened JSON → `[meta, file]`; blank/optional path → `[]` |
| Synchronized dispatch + errors-as-result | `ModuleToolBridge.call` | LLM can recover from bad args without aborting the run |
| Read outputs, skip `eval`/`topic`, serialize | `ModuleToolBridge.callSpec` + `isTopicOutput`/`isEvalOutput` | topic-source channels never bind → would hang |
| Inline small structured outputs vs path handles | `ModuleToolBridge.serializeValue` → `ToolOutputReader.readOrHandle` | enables LLM reasoning (gate) vs chaining (handles) |
| Per-tool `moduleDir`/binding isolation | `AgentDef.compileModuleProcess` + `ScriptLoaderV2.setMainBinding` | otherwise all tools share one container |
| Tool descriptor/schema from registry | `ModuleMetadataToolSchema` + `JsonSchemaMapper` | high-level instructions; self-describing tools |
| The LLM tool-calling loop | `LangChainAgentRunner.runWithTools` | advertise tools → dispatch → feed result → repeat |
| Poison-on-complete | `AgentDef.applyAgentOperator` (`afterStop`) → `ModuleToolBridge.close` | clean dataflow termination |

---

## 9. Bugs fixed in this work

1. **Inline small structured outputs** (`ToolOutputReader` + `agent.maxToolOutputInlineSize`) —
   surface small `.json`/`.tsv` tool outputs to the LLM; bulk data stays a handle.
2. **Per-tool `moduleDir` isolation** — each tool module gets its own `ScriptBinding`
   (else all tools share one container).
3. **Skip topic-routed `versions` outputs via the ProcessDef** — the meta.yml `type` is
   unreliable; reading a topic-source `.val` hangs.
4. **V2 guard** — `isTopicOutput` guards `BaseOutParam` (typed `ProcessOutput.getChannelTopicName()`
   throws).
5. **Blank/optional file inputs** — a `""` path arg → not provided (`[]`), never `file("")`.
6. **Capture output read channels at wiring time** — fixes the multi-output dispatch deadlock.

## 10. Known limitation (open)

A **genuinely absent optional output** (a tool that produces *no* file for an `optional: true`
output) is still unhandled: `TaskProcessor.bindOutputsV1` skips binding it, so the channel never
gets a value for that call and the per-call `.val` blocks until process termination (which the
blocked dispatch prevents). The current example tools always produce their declared outputs, so
this does not arise there. A robust fix needs either a non-terminating "empty for this task"
sentinel on the output channel in core, or a bridge-local guard — a core dataflow-semantics
decision deferred pending direction.
