#!/usr/bin/env node

import process from "node:process";
import {
  createAgentSession,
  DefaultResourceLoader,
  defineTool,
  ModelRuntime,
  SessionManager,
} from "@earendil-works/pi-coding-agent";

// Bumped to 2 when `spec.nativeToolNames` was added. A new field is only backward-compatible
// when IGNORING it is safe, and this one is not: a v1 harness reads no native names and finds no
// `filesystem` descriptor either (the driver stopped sending one), so it would run an agent whose
// instruction names file and shell tools it was never given, and fail as a confusing model error
// rather than a version mismatch. Refusing the frame turns that into one legible line.
const PROTOCOL_VERSION = 2;
const pendingToolCalls = new Map();
let activeInvocation;
let session;
let terminal = false;
let stdinBuffer = "";

function send(message) {
  process.stdout.write(`${JSON.stringify(message)}\n`);
}

// stdout is the exclusive JSONL protocol channel: only send() may write to it.
// Redirect every diagnostic console channel to stderr so a stray log line from
// the Pi SDK or a transitive dependency cannot corrupt or inject a frame.
console.log = (...args) => process.stderr.write(`${args.join(" ")}\n`);
console.info = console.log;
console.debug = console.log;
console.warn = console.log;
// console.error already targets stderr; leave it untouched.

function fail(error, code = "runner_error") {
  if (terminal) return;
  terminal = true;
  send({
    type: "error",
    invocationId: activeInvocation,
    code,
    message: error instanceof Error ? error.message : String(error),
  });
}

function composeSystemPrompt(spec) {
  const parts = [];
  if (spec.instruction) parts.push(spec.instruction);
  if (spec.goal) parts.push(`Goal:\n${spec.goal}`);
  if (spec.skills?.length) {
    parts.push(
      "Available skills:\n" +
        spec.skills
          .map((skill) => `- ${skill.name}: ${skill.description}`)
          .join("\n") +
        "\n\nUse activate_skill when a skill is relevant, then follow its instructions. " +
        "Use read_skill_resource for any bundled resource named by the activated skill.",
    );
  }
  if (spec.outputSchema) {
    parts.push(
      "You MUST finish by calling final_answer exactly once with arguments matching its schema. " +
        "Do not finish with ordinary assistant text.",
    );
  }
  return parts.join("\n\n");
}

function composeUserPrompt(spec) {
  if (!spec.inputJson) return spec.prompt ?? "";
  return `${spec.prompt ?? ""}\n\nInput:\n${spec.inputJson}`;
}

// The task work dir the driver assigned, or this process' own -- the resource loader, the session
// and its manager must all agree on it, so it is spelled once.
function workDirOf(spec) {
  return spec.workDir || process.cwd();
}

function brokerTool(name, toolCallId, params, signal) {
  return new Promise((resolve, reject) => {
    if (signal?.aborted) {
      reject(new Error(`Tool ${name} was cancelled`));
      return;
    }
    const onAbort = () => {
      pendingToolCalls.delete(toolCallId);
      reject(new Error(`Tool ${name} was cancelled`));
    };
    signal?.addEventListener("abort", onAbort, { once: true });
    pendingToolCalls.set(toolCallId, {
      resolve: (value) => {
        signal?.removeEventListener("abort", onAbort);
        resolve(value);
      },
      reject: (error) => {
        signal?.removeEventListener("abort", onAbort);
        reject(error);
      },
    });
    send({
      type: "tool_call",
      invocationId: activeInvocation,
      callId: toolCallId,
      name,
      arguments: params,
    });
  });
}

// The `maxIterations` budget, applied to every tool this harness defines EXCEPT `final_answer`:
// that one neither counts nor checks, so an agent that reaches the cap exactly can still return
// its answer rather than failing on the way out.
//
// Two details are load-bearing and invisible to the scripted stub:
//   - the delegation forwards the FULL argument list, because the wrapped `execute`s do not share
//     a signature -- a brokered tool takes `(toolCallId, params, signal)` and passes the abort
//     signal on, while the skill tools take `(_toolCallId, params)`. Dropping `signal` would leave
//     a cancelled invocation waiting out its RPC timeout instead of aborting;
//   - the counter is incremented BEFORE the comparison, which is `>`. Checking first, or using
//     `>=`, shortens every agent's tool budget by one call.
function cappedTool(spec, state, definition) {
  return {
    ...definition,
    async execute(...args) {
      state.toolTurns += 1;
      if (state.toolTurns > spec.maxIterations) {
        throw new Error(
          `Agent exceeded the maximum number of tool-call iterations (${spec.maxIterations})`,
        );
      }
      return definition.execute(...args);
    },
  };
}

// Every descriptor in `toolSpecs` is BROKERED, unconditionally: the driver owns the tools it
// describes, and the runner's own tools never arrive as descriptors -- they arrive as bare names
// in `spec.nativeToolNames` and are enabled from the SDK's builtins (see the session allowlist
// below). There is deliberately no per-descriptor locality flag: a mis-set one would execute a
// container-side tool in the driver JVM, and it would oblige this harness to reimplement what the
// SDK already ships.
function nextflowTools(spec, state) {
  const result = [];
  for (const descriptor of spec.toolSpecs ?? []) {
    result.push(
      defineTool(
        cappedTool(spec, state, {
          name: descriptor.name,
          label: descriptor.name,
          description: descriptor.description || `Run Nextflow tool ${descriptor.name}`,
          promptSnippet: descriptor.description || `Run Nextflow tool ${descriptor.name}`,
          parameters: descriptor.inputSchema,
          executionMode: "sequential",
          async execute(toolCallId, params, signal) {
            const text = await brokerTool(descriptor.name, toolCallId, params, signal);
            return { content: [{ type: "text", text }], details: {} };
          },
        }),
      ),
    );
  }
  return result;
}

function skillTools(spec, state) {
  if (!spec.skills?.length) return [];
  const byName = new Map(spec.skills.map((skill) => [skill.name, skill]));
  return [
    defineTool(
      cappedTool(spec, state, {
        name: "activate_skill",
        label: "Activate skill",
        description: "Load the complete instructions for an available skill.",
        promptSnippet: "Load the instructions for a relevant available skill",
        parameters: {
          type: "object",
          properties: { name: { type: "string" } },
          required: ["name"],
          additionalProperties: false,
        },
        executionMode: "sequential",
        async execute(_toolCallId, params) {
          const skill = byName.get(params.name);
          if (!skill) throw new Error(`Unknown skill: ${params.name}`);
          state.activeSkills.add(skill.name);
          return { content: [{ type: "text", text: skill.content }], details: {} };
        },
      }),
    ),
    defineTool(
      cappedTool(spec, state, {
        name: "read_skill_resource",
        label: "Read skill resource",
        description: "Read a bundled resource from an activated skill.",
        promptSnippet: "Read a bundled resource from an activated skill",
        parameters: {
          type: "object",
          properties: {
            skill: { type: "string" },
            path: { type: "string" },
          },
          required: ["skill", "path"],
          additionalProperties: false,
        },
        executionMode: "sequential",
        async execute(_toolCallId, params) {
          if (!state.activeSkills.has(params.skill))
            throw new Error(`Skill is not active: ${params.skill}`);
          const skill = byName.get(params.skill);
          const resource = skill?.resources?.find((item) => item.relativePath === params.path);
          if (!resource) throw new Error(`Unknown skill resource: ${params.skill}/${params.path}`);
          return { content: [{ type: "text", text: resource.content }], details: {} };
        },
      }),
    ),
  ];
}

function structuredOutputTool(spec, state) {
  if (!spec.outputSchema) return [];
  return [
    defineTool({
      name: "final_answer",
      label: "Final answer",
      description: "Return the final structured answer and terminate the agent.",
      promptSnippet: "Return the final structured answer and terminate",
      promptGuidelines: ["Use final_answer exactly once as the last action for structured output."],
      parameters: spec.outputSchema,
      executionMode: "sequential",
      async execute(_toolCallId, params) {
        state.finalOutput = JSON.stringify(params);
        return {
          content: [{ type: "text", text: "Structured answer accepted" }],
          details: {},
          terminate: true,
        };
      },
    }),
  ];
}

async function run(start) {
  const spec = start.spec;
  activeInvocation = start.invocationId;
  if (start.protocolVersion !== PROTOCOL_VERSION)
    throw new Error(`Unsupported protocol version: ${start.protocolVersion}`);
  if (!spec?.model) throw new Error("Agent `model` directive is required");

  const slash = spec.model.indexOf("/");
  if (slash <= 0 || slash === spec.model.length - 1)
    throw new Error(`Invalid model identifier: ${spec.model}; expected provider/model`);
  const provider = spec.model.slice(0, slash);
  const modelId = spec.model.slice(slash + 1);
  const modelRuntime = await ModelRuntime.create({ allowModelNetwork: false });
  // A resolved endpoint RETARGETS the provider catalog: registerProvider with no `models`
  // rewrites baseUrl on every built-in model of the provider, and is the only seam the SDK
  // sanctions. The endpoint must NOT be spread onto the resolved model object instead: the
  // catalog owns Model instances and re-resolves them by id (agent-session refreshes
  // `state.model` from the registry, session restore does the same), so an ad-hoc field on a
  // model is fragile. Registering an unknown provider id is harmless -- it composes an empty
  // catalog, so the model lookup below still reports `Unknown Pi model`.
  // The retargeted endpoint must speak the wire protocol of the provider being retargeted: the
  // built-in `openai` provider is wired to the OpenAI RESPONSES api, so it posts /v1/responses.
  // An OpenAI mirror or a Responses-capable gateway works; a chat/completions-only server
  // (Ollama, llama.cpp, default vLLM) does not, and is out of scope here.
  if (spec.baseUrl) modelRuntime.registerProvider(provider, { baseUrl: spec.baseUrl });
  const model = modelRuntime.getModel(provider, modelId);
  // Only catalog model ids resolve: `baseUrl` retargets known ids at a compatible endpoint,
  // it does not add models. Registering an arbitrary id needs a `models: [...]` entry that
  // REPLACES the provider catalog and mandates per-model metadata, so it stays out of scope
  // and this failure is the documented boundary. Checked BEFORE the credential is installed
  // so an unknown id never surfaces as a confusing auth error.
  if (!model) throw new Error(`Unknown Pi model: ${spec.model}`);
  // The credential is delivered OUT OF BAND (beside `spec`, never inside it) and installed as
  // a pi runtime API key: it is held in memory for this process only, so it neither enters
  // this process' environment nor pi's on-disk auth store.
  // A runtime key OWNS the provider -- runtime-credentials.js returns the override before the
  // store, and auth/resolve.js consults the ambient environment only when nothing is stored --
  // so the driver sends a key here ONLY when it is a REAL credential it resolved in this
  // provider's own namespace (`agent.apiKey`/`NXF_AGENT_API_KEY`, or `<PROVIDER>_API_KEY` for an
  // endpoint that belongs to that provider). Never its no-credential placeholder: that one owns
  // nothing, and installing it would shadow whatever the container was given out of band. When
  // nothing arrives -- no credential resolved, `agent.rpc.tls = false`, or a key the driver's
  // endpoint gate withheld (all three make the driver warn rather than fail, precisely because
  // of what follows) -- pi's own resolution (its store, then provider variables such as
  // ANTHROPIC_API_KEY delivered by the `env` scope or a Kubernetes Secret) is left untouched.
  if (start.apiKey) await modelRuntime.setRuntimeApiKey(provider, start.apiKey);

  const state = { toolTurns: 0, activeSkills: new Set(), finalOutput: undefined };
  const customTools = [
    ...nextflowTools(spec, state),
    ...skillTools(spec, state),
    ...structuredOutputTool(spec, state),
  ];
  const systemPrompt = composeSystemPrompt(spec);
  const loader = new DefaultResourceLoader({
    cwd: workDirOf(spec),
    agentDir: workDirOf(spec),
    noExtensions: true,
    noSkills: true,
    noPromptTemplates: true,
    noThemes: true,
    noContextFiles: true,
    systemPromptOverride: () => systemPrompt,
    appendSystemPromptOverride: () => [],
  });
  await loader.reload();

  // The allowlist is the WHOLE tool gate, and it is exhaustive by construction: the tools this
  // harness defines (brokered Nextflow tools, skills, final_answer) plus the SDK builtins the
  // agent selected through the `fs:`/`shell:` families. Anything absent is disabled -- an
  // allowlist is the only thing pi consults when one is given (`allowedToolNames = options.tools`
  // in dist/core/sdk.js), which is also why `noTools` is gone: `noTools: "builtin"` was dead
  // whenever `tools` was set, and it now reads as the opposite of what this line does, since
  // enabling the resolved builtins is the entire point of the runner-native half of the split.
  // The builtins are built for the session `cwd`, i.e. the task work dir passed below, and the
  // container is their outer bound.
  // NOTE a builtin is executed by the SDK, so it does NOT pass through `state.toolTurns`: the
  // `maxIterations` budget bounds the tools defined here (brokered, skills, final_answer) and not
  // a `read`/`grep` loop inside the container. Bounding those needs a counter driven off the
  // session event stream, which is a separate change.
  const created = await createAgentSession({
    cwd: workDirOf(spec),
    modelRuntime,
    model,
    tools: [...customTools.map((tool) => tool.name), ...(spec.nativeToolNames ?? [])],
    customTools,
    resourceLoader: loader,
    sessionManager: SessionManager.inMemory(workDirOf(spec)),
  });
  session = created.session;

  let currentText = "";
  let lastAssistantText = "";
  // Last provider/transport failure seen on the event stream. The SDK reports a failed
  // model call as an assistant message with `stopReason: "error"`, an `errorMessage`, and
  // EMPTY content, then auto-retries. Without capturing it, an exhausted retry chain
  // surfaces only as "no output", which misattributes a provider outage to the model
  // declining to answer and hides the request id needed to report it.
  let lastProviderError = "";
  session.subscribe((event) => {
    if (event.message?.errorMessage) lastProviderError = event.message.errorMessage;
    if (event.errorMessage) lastProviderError = event.errorMessage;
    if (event.finalError) lastProviderError = String(event.finalError.message ?? event.finalError);
    if (event.type === "message_start" && event.message?.role === "assistant") currentText = "";
    if (event.type === "message_update") {
      const update = event.assistantMessageEvent;
      if (update?.type === "text_delta") currentText += update.delta;
      if (update?.type === "thinking_delta")
        send({ type: "trace", invocationId: activeInvocation, event: "thinking", text: update.delta });
    }
    if (event.type === "message_end" && event.message?.role === "assistant")
      lastAssistantText = currentText;
    if (event.type === "tool_execution_start")
      send({
        type: "trace",
        invocationId: activeInvocation,
        event: "tool_start",
        name: event.toolName,
        callId: event.toolCallId,
      });
    if (event.type === "tool_execution_end")
      send({
        type: "trace",
        invocationId: activeInvocation,
        event: "tool_end",
        name: event.toolName,
        callId: event.toolCallId,
        isError: event.isError,
      });
  });

  try {
    await session.prompt(composeUserPrompt(spec));
    if (spec.outputSchema && state.finalOutput == null) {
      await session.prompt(
        "Your previous response did not satisfy the required machine-readable output contract. " +
          "Call final_answer now with the complete answer. Do not respond with ordinary text.",
      );
    }
    const output = spec.outputSchema ? state.finalOutput : lastAssistantText;
    if (output == null || output.trim() === "") {
      const reason = spec.outputSchema ? "Pi did not call final_answer" : "Pi returned no final assistant text";
      throw new Error(lastProviderError ? `${reason}; last provider error: ${lastProviderError}` : reason);
    }
    terminal = true;
    send({
      type: "complete",
      invocationId: activeInvocation,
      output,
      resolvedModel: `${model.provider}/${model.id}`,
    });
  } finally {
    session.dispose();
    session = undefined;
  }
}

async function handle(message) {
  if (message.type === "start") {
    if (activeInvocation) throw new Error("Runner accepts exactly one invocation");
    await run(message);
    return;
  }
  if (message.type === "tool_result") {
    if (message.invocationId !== activeInvocation) throw new Error("Mismatched invocationId");
    const pending = pendingToolCalls.get(message.callId);
    if (!pending) throw new Error(`Unknown tool call result: ${message.callId}`);
    pendingToolCalls.delete(message.callId);
    if (message.isError) pending.reject(new Error(message.result));
    else pending.resolve(message.result);
    return;
  }
  if (message.type === "cancel") {
    session?.abort();
    throw new Error(message.reason || "Agent invocation cancelled");
  }
  throw new Error(`Unknown host message type: ${message.type}`);
}

process.stdin.setEncoding("utf8");
process.stdin.on("data", (chunk) => {
  stdinBuffer += chunk;
  for (;;) {
    const newline = stdinBuffer.indexOf("\n");
    if (newline < 0) break;
    const line = stdinBuffer.slice(0, newline).replace(/\r$/, "");
    stdinBuffer = stdinBuffer.slice(newline + 1);
    if (!line) continue;
    let message;
    try {
      message = JSON.parse(line);
    } catch (error) {
      fail(error, "invalid_json");
      continue;
    }
    handle(message).catch((error) => fail(error));
  }
});
process.stdin.on("end", () => {
  if (!terminal) fail(new Error("Host closed protocol input"), "unexpected_eof");
});
process.on("uncaughtException", (error) => fail(error));
process.on("unhandledRejection", (error) => fail(error));

send({ type: "ready", protocolVersion: PROTOCOL_VERSION });
