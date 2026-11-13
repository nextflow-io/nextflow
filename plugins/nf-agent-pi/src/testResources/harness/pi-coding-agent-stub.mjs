// Copyright 2013-2026, Seqera Labs
// SPDX-License-Identifier: Apache-2.0
//
// Test double for `@earendil-works/pi-coding-agent`, the only package harness/runner.mjs
// imports. PiHarnessProtocolTest installs it into a throwaway `node_modules` beside a copy of
// the harness, so the REAL harness can be driven end to end with no `npm ci`, no network and
// no provider credentials -- the runtime now lives in the container image, and the Gradle
// build must not need npm to test it.
//
// The scripted model behaviour comes from PI_STUB_PLAN: a JSON array of TURNS, one consumed
// per `session.prompt()` (a second turn is the harness's corrective re-prompt). Each turn is
// an array of steps:
//
//   {"text": "..."}                     emit assistant text
//   {"tool": "<name>", "args": {...}}   call a tool the harness registered
//   {"tool": ..., "echo": true}         ... and make its result the assistant's answer
//   {"providerError": "..."}            a failed model call: an assistant message carrying an
//                                       errorMessage and NO content, as the SDK reports it
//   {"runtimeState": true}              answer with the ModelRuntime calls the harness made and
//                                       the session options it built, so a spec can assert on
//                                       provider retargeting, credentials and the tool allowlist

import process from "node:process";

const turns = JSON.parse(process.env.PI_STUB_PLAN ?? "[]");

/**
 * What the harness asked of the ModelRuntime before the session started, plus the tool selection
 * it handed to createAgentSession. Recorded rather than acted on: these are how an endpoint, a
 * credential and the enabled tool set reach the model, and none of them is observable from the
 * protocol frames alone.
 *
 * `session.tools` is the SDK's ALLOWLIST -- the runner-native half of the agent's tools is enabled
 * by naming builtins in it, so it is the only place a spec can see that half at all (a builtin has
 * no descriptor and never appears as a customTool).
 */
const runtimeState = { providers: {}, apiKeys: {}, session: {} };

/** The SDK returns the definition enriched; the harness only ever reads `.name` off it. */
export function defineTool(definition) {
  return definition;
}

export class ModelRuntime {
  static async create() {
    return new ModelRuntime();
  }
  /** Retargets the catalog of a provider; with no `models` it only rewrites the endpoint. */
  registerProvider(provider, options) {
    runtimeState.providers[provider] = options;
  }
  /** An in-memory credential that OWNS the provider, taking precedence over the auth store. */
  async setRuntimeApiKey(provider, apiKey) {
    runtimeState.apiKeys[provider] = apiKey;
  }
  getModel(provider, id) {
    return id === "no-such-model" ? null : { provider, id };
  }
}

export class DefaultResourceLoader {
  constructor(options) {
    this.options = options;
  }
  async reload() {
    // exercise the overrides the harness passes, so a rename is not silently ignored
    this.systemPrompt = this.options.systemPromptOverride();
    this.appended = this.options.appendSystemPromptOverride();
  }
}

export const SessionManager = {
  inMemory: (cwd) => ({ cwd }),
};

export async function createAgentSession(options) {
  const { customTools } = options;
  // The harness redirects console.log to stderr before it gets here, so a chatty SDK cannot
  // corrupt the JSONL protocol channel. If that redirect ever goes away this line lands on
  // stdout and every spec fails to parse a frame -- which is the point.
  console.log("pi stub: this console.log must not reach stdout");

  runtimeState.session = {
    cwd: options.cwd,
    tools: options.tools,
    // recorded so a spec can assert on its ABSENCE: the allowlist wins over it in the real SDK
    noTools: options.noTools ?? null,
    customTools: customTools.map((tool) => tool.name),
  };
  const tools = new Map(customTools.map((tool) => [tool.name, tool]));
  const listeners = [];
  const emit = (event) => listeners.forEach((listener) => listener(event));
  let callSeq = 0;
  let aborted = false;

  const say = (text) => {
    emit({ type: "message_start", message: { role: "assistant" } });
    emit({ type: "message_update", assistantMessageEvent: { type: "text_delta", delta: text } });
    emit({ type: "message_end", message: { role: "assistant" } });
  };

  const session = {
    subscribe(listener) {
      listeners.push(listener);
    },
    abort() {
      aborted = true;
    },
    dispose() {},
    async prompt() {
      for (const step of turns.shift() ?? []) {
        if (aborted) return;
        if (step.text !== undefined) {
          say(step.text);
          continue;
        }
        if (step.runtimeState !== undefined) {
          say(JSON.stringify(runtimeState));
          continue;
        }
        if (step.providerError !== undefined) {
          emit({
            type: "message_end",
            message: { role: "assistant", errorMessage: step.providerError },
          });
          continue;
        }
        const tool = tools.get(step.tool);
        if (!tool) throw new Error(`pi stub: the harness registered no tool named ${step.tool}`);
        const callId = `call-${++callSeq}`;
        emit({ type: "tool_execution_start", toolName: tool.name, toolCallId: callId });
        const result = await tool.execute(callId, step.args ?? {}, undefined);
        emit({
          type: "tool_execution_end",
          toolName: tool.name,
          toolCallId: callId,
          isError: false,
        });
        if (step.echo) say(result.content.map((part) => part.text).join(""));
        if (result.terminate) return;
      }
    },
  };

  return { session };
}
