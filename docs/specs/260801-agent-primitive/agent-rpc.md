# Agent RPC — canonical task execution and the driver broker

- Status: implemented for direct driver connections; parts of the production hardening remain
  proposed and are marked as such
- Scope: lowering an agent to a normal executor-submitted task, the protocol between that task and the
  driver, transport security, driver-address inference, and scalability
- Companion documents: [agent design](agent-design.md) (the primitive),
  [agent runners](agent-runners.md) (which runners use this transport)

## 1. The problem this solves

The first out-of-process runner still ran its loop from inside the Nextflow JVM: the agent was a
`TaskRun`, but its body was a Groovy closure that spawned a child process. Remote executors submit
`SCRIPTLET` tasks, so a Groovy body is always executed locally — the harness could be put in a local
container, but it could never be submitted by Kubernetes, Slurm or Batch.

The obvious alternative, **checkpointing at every tool boundary**, makes each segment independently
schedulable at the cost of a pod start, a runtime initialization and a session restore *per turn*.
Convergence loops make that overhead impossible to ignore.

So the target is:

```
one logical agent invocation
=  one canonical local or remote Nextflow task
+  one proxy holding one long-lived outbound RPC stream
+  one harness child process on JSONL stdio
+  zero or more canonical Nextflow tool tasks
```

```
Nextflow driver                        Canonical local/remote agent task
┌──────────────────────┐              ┌───────────────────────────┐
│ JVM Agent RPC broker │◄══ gRPC ════►│ Go agent-rpc task command │
│ Tool dispatcher      │              │            ↕ JSONL stdio  │
│ Workflow runtime     │              │ Pi or third-party harness │
└──────────┬───────────┘              └───────────────────────────┘
           ├──► SKESA task on Kubernetes
           └──► other canonical tasks
```

The feature is **not** a new agent-specific executor. It is lowering an agent invocation to a normal
`SCRIPTLET` `TaskRun` submitted through the existing `ExecutorFactory` / `Executor` / `TaskHandler` /
staging / trace / retry / cancellation lifecycle.

## 2. Design principles

1. The agent task is submitted and monitored through the **unmodified** executor lifecycle.
2. The task-side proxy initiates all network connections; workers need no inbound reachability.
3. RPC carries intermediate control messages, **not task completion authority** — the `TaskHandler`
   remains the authoritative source of task status, and a disconnected stream is not by itself
   evidence that a pod failed.
4. `ToolDispatcher` and workflow runtime objects never leave the JVM.
5. Module tools execute as real Nextflow tasks using **their own** executor and container
   configuration.
6. Agent placement and tool placement are independent decisions. All four combinations of
   local/remote agent × local/remote tool are architecturally supported, and three were validated.
7. Tool delivery is at least once; execution is idempotent by invocation, attempt and call ID.
8. No runner may invoke a tool outside the allowlist issued for its task.
9. The Nextflow installation does not depend on a platform-specific binary, and a harness needs
   neither Java nor gRPC.

Agent tasks resolve placement and resources from the `agent` scope, defaulting to
`agent.executor = 'local'` — never from `process` defaults or selectors. Wave, Fusion,
container-engine and work-directory scopes remain shared session services. Kubernetes `pod`
directives are deferred: their additive/list semantics need a dedicated configuration model rather
than a scalar pass-through.

## 3. Transport

```
Canonical task ↔ driver:      gRPC bidirectional streaming over TLS
Proxy ↔ harness:              JSONL over stdio
Unit-test compatibility path: JSONL over stdio, no proxy
```

```protobuf
service AgentBroker {
  rpc Connect(stream RunnerMessage) returns (stream HostMessage);
}
```

The semantic protocol is **independent of gRPC**. Invocation identity, call idempotency, tool
authorization and terminal-state validation are shared by every transport; JVM-side adapters feed one
invocation registry, and the proxy translates frames without owning workflow semantics. JSON-RPC over
WebSocket may later be added as a compatibility adapter for third-party harnesses, but only against
the same state machine — not with different delivery, authorization or recovery semantics.

Frames currently in use:

| Direction | Messages |
|---|---|
| harness → driver | `connect`, `tool_call`, `trace`, `complete`, `error` |
| driver → harness | `ready`, `start` (the portable invocation spec, including the credential), `tool_result`, `cancel` |

Every message after `connect` carries the `invocationId`; tool frames carry a unique `callId`. The
driver **memoizes completed `callId` results** for the invocation's lifetime, so a duplicated protocol
frame cannot execute an expensive scientific module twice. Unknown message types, malformed JSON,
oversized lines, reused call IDs with changed payloads and mismatched invocation IDs are fatal
protocol errors. Stderr is diagnostic output and is never parsed as protocol data.

**Deferred production vocabulary.** Durable identity (`workflowRunId`, `taskAttempt`, `sequence`),
`heartbeat`, `ack`, and an acknowledged terminal state in which the proxy writes declared result and
manifest artifacts before exiting — so that completion behaviour is uniform across harness
implementations — are designed but not implemented. Today the proxy prints the `complete` payload on
stdout and the driver decodes the last non-empty line.

## 4. Lifecycle

Before submission the driver assigns an invocation ID, registers it with the broker, issues a
short-lived capability token, retains the portable agent specification driver-side, and submits the
task through the executor the `agent` scope selects. The generated script is equivalent to:

```bash
/usr/local/bin/agent-rpc \
    --endpoint <host>:<port> \
    --invocation <INVOCATION_ID> \
    --token <ONE_USE_TOKEN> \
    --fingerprint <CERT_SHA256> \
    -- node /opt/nf-agent-pi/runner.mjs
```

```
harness ↔ JSONL ↔ proxy                gRPC/TLS          driver broker
   ├── ready
   │                     connect ──────────────►   (invocationId + one-use token)
   │◄──────────────────── start/spec ─────────┤
   ├── tool_call (callId, name, args) ───────►├── schedule a Nextflow task
   │◄──────────────────── tool_result ────────┤
   ├── trace ────────────────────────────────►│
   ├── complete ─────────────────────────────►│
   └── print the complete payload, exit 0
```

## 5. Authorization and transport security

### 5.1 What the capability is

32 bytes of `SecureRandom`, base64url-encoded, **single use** — consumed atomically on connect, so an
observed token cannot be replayed and a second stream cannot race a connected one. Every subsequent
frame is authorized by invocation identity instead. The token is never logged, and provider keys are
redacted from diagnostics. Dispatch enforces a per-invocation tool allowlist, call-ID deduplication
and an atomic call ceiling.

The authorization model was never the problem. **Lifetime** and **transport** were.

### 5.2 Capability lifetime

The capability is registered when the task **script is generated**, which `TaskProcessor` does long
before it submits. The original budget of `max(requestTimeout, 30) + 60` therefore started ticking
during queueing, and any executor that queues — Kubernetes under pressure, Slurm, AWS Batch — would
routinely blow through 180 seconds and fail with `Invalid agent RPC invocation identity or token`: a
security-shaped message for what is really scheduling latency, pointing the operator in entirely the
wrong direction.

Fix: `agent.rpc.capabilityTimeout`, default **1 hour**, governing the pre-connect wait only. Because
the value is a `Duration`, a bare number is milliseconds — `= 3600` is three seconds, not an hour —
which the user guide calls out in two places.

Two rejected alternatives, both worth recording:

- **Revoking the capability when its task reaches a terminal state.** Theoretically tidy, but it makes
  a plugin-hosted broker observe `TaskRun` lifecycle across a module boundary and fails badly in both
  directions: a missed hook leaks capabilities for the run, a misfiring one revokes mid-stream and
  produces exactly the spurious security error this fix exists to remove. A constant cannot misfire.
- **A post-connect timer.** The design originally kept the old budget and re-armed it at connect. That
  budget is sized for *one* model request, and a connected agent legitimately makes many — with
  `maxIterations` 20 and `requestTimeout` 120 s, ~40 minutes of model time alone, plus every
  `tool_call` blocking on a real task that queues behind the pipeline's own work. A 180 s post-connect
  deadline would kill essentially every non-trivial agent, reintroducing the same error at the other
  end of the lifecycle; re-arming at a larger constant only moves the guess. **Nothing is armed
  post-connect.** What such a timer actually covered — a peer that vanishes without a FIN, an
  OOM-killed pod, a reclaimed spot instance — is covered by gRPC keepalive, which notices a dead node
  in ~80 s without capping a live one.

`ScheduledThreadPoolExecutor` defaults `removeOnCancelPolicy` to `false`, so the pool is built
explicitly with it enabled — otherwise a *successful* invocation's cancelled deadline stays in the
delay queue, pinning its whole `AgentRunnerRequest` for an hour after it finished. A `-resume` cache
hit registers a capability nothing ever consumes (registration happens before the cache is consulted),
so unconnected capabilities are counted and reported once at shutdown; without that the retention
would be invisible rather than absent.

### 5.3 TLS with a pinned ephemeral certificate

The broker generates an EC P-256 key pair and a self-signed certificate **in memory, once per run**,
serves gRPC over TLS with it, and passes the certificate's SHA-256 fingerprint to the task, which
connects only on an exact match. This is the SSH-known-hosts pattern, and the key move is that **a
fingerprint is a public commitment, not a secret** — putting it in `.command.sh` is harmless.

With no PKI, no CA, no certificate files and nothing to rotate, it buys:

- **payload confidentiality** — the start frame carries model, instruction, goal, prompt, input JSON,
  output schema, tool and skill specs and the credential, and every tool argument and result crosses
  the same link. For this domain that is sample identifiers, phenotype text, file paths and
  intermediate results;
- **server authentication** — a spoofed endpoint cannot present a matching certificate, closing driver
  impersonation. The inbound direction was already defended by the unguessable single-use token; the
  *outbound* trust direction had nothing;
- defence in depth for the token, which no longer crosses the network in cleartext.

BouncyCastle supplies the certificate object — the JDK has no public API for generating a self-signed
X.509 — at the version `nf-k8s` already pins, declared independently because PF4J isolates plugin
classloaders. It is **not** registered as a JVM-wide security provider. Netty's bundled
`SelfSignedCertificate` was rejected for living at a shaded package path that is not stable across
gRPC upgrades; shelling out to `keytool` for being a process spawn per run. The ~10 MB this adds back
to a plugin the packaging change had just taken from 186 MB to under 1 MB is an ordinary pinned Maven
dependency — pure bytecode, every architecture, no extraction, no toolchain — so none of the four
defects that work fixed returns; only the byte count partly rebounds, and bytes were the symptom.
Slimming it is not an option: BouncyCastle jars are signed, and shading breaks the signature.

`agent.rpc.tls = false` remains as a debugging escape hatch, with a warning, and withholds the
credential.

### 5.4 The accepted residual: the token on argv

The token is on the task command line, so it is written into `.command.sh` and persists in the work
directory — commonly group-readable on a shared filesystem, and subject to bucket ACLs on object
storage.

Calibrate it honestly. Anyone who can read `.command.sh` can already read the staged inputs,
`.command.env` and the outputs — the pipeline's data is already theirs. The token's *marginal* harm is
live invocation hijack, driver-side tool dispatch within the invocation's allowlist, and — since the
driver answers a connection with the start frame — **disclosure of the provider API key**. That last
one is a consequence of in-band credential delivery that the credential design did not consider.

Certificate pinning does not help here: it is a client-side check, so a thief simply connects to the
real driver and validates. **The lifetime fix is the load-bearing mitigation** — a capability that is
single-use and dies with its window means the durable copy is almost always already dead.

What *was* fixed is the token escaping the work directory. `TaskHandler.getTraceRecord` set
`record.script` verbatim, and a `TraceRecord` is persisted in the resume cache and POSTed to Seqera
Platform — turning a work-dir-local secret into a stored and transmitted one, for the life of the
cache. `TaskRun.getTraceScript()` now redacts `--token` **for agent tasks only**, guarded by the same
check lineage already uses to omit the agent script entirely. Every non-agent task gets a
byte-identical value, and the executed `.command.sh` keeps the real token.

Rejected for token delivery: **Nextflow secrets** (writing per-invocation machine state into a
user-facing secret store, racy and invasive); an **environment variable** (Nextflow writes task env
into `.command.run`, the same durable artifact with extra steps); and **mTLS with a per-invocation
client certificate** (strictly better cryptographically, but the private key has exactly the delivery
problem the token has). A `--token-file` with mode 0600 is a contained additive follow-up if a
security review requires it — the gain is POSIX-only and **zero** on object-storage work dirs, which
is why it was not shipped.

`agent.rpc.bindAddress` was specified and then dropped. Its impersonation half is closed by pinning
outright; its exposure half is not what a bind address fixes, since practically every agent task now
reaches the driver from *outside* the driver's network namespace. What would remain is an option whose
wrong value fails silently in the wrong direction. The server binds every interface, and the code says
so explicitly rather than leaving it to be inferred.

## 6. Reaching the driver

Because the runner is container-only, `agent.rpc.remoteHost` became load-bearing on every run — and
its old default, `host.docker.internal`, resolves only on Docker Desktop. Advertising an address
nothing can reach surfaces as a connect timeout *inside the task*, an hour later, possibly after paying
for model calls.

The address is now **inferred** from what the driver already knows, and the cases that genuinely
cannot work are rejected **before ignition**, each with its own message.

| # | Condition | Address |
|---|---|---|
| R1 | `agent.rpc.remoteHost` set | that value |
| R2 | `NXF_AGENT_RPC_REMOTE_HOST` set | that value |
| R3 | local executor, `docker` or `podman` | `host.docker.internal` / `host.containers.internal` |
| R4 | local executor, `apptainer` or `singularity` | `127.0.0.1` |
| R5 | local executor, `smolvm` (networking enabled) or `apple-container` | the driver's default-route address |
| R6 | everything else — any other engine on the driver host, grid, Kubernetes, cloud batch | the driver's default-route address |

R4 is a property of the engine, not a guess: apptainer and singularity create no network namespace,
so the driver's loopback *is* the container's. An explicit `--network host` does not create a rung of
its own. Every other engine on the driver host — Shifter, Charliecloud, Sarus, anything unrecognized
— falls through to R6 rather than to loopback: an unknown engine still containerizes, it simply names
no address for the driver host.

Only three configurations are rejected: **no container engine enabled at all** (nothing is
containerized, so there is no task to reach back), **`smolvm.network = false`** (the microVM is
created with no network), and **no usable address** (the default route yields loopback, a wildcard,
a link-local address or nothing, and the interfaces do not settle on one unambiguous answer).

Two measured facts drive this. **The guest's own default gateway is never a candidate** — for Docker
it is the `docker0` bridge inside Docker Desktop's Linux VM, and dialling it fails. And **the host's
default-route source address covers every engine that could be tested**, including both VM-isolated
ones, so a single advertised address suffices.

The outbound address is the local address the kernel picks for the default route
(`DatagramSocket.connect(…)` then `getLocalAddress()`) — a routing-table lookup, no packet sent. It is
memoized per session, since it is a property of the host rather than of the agent. Loopback, wildcard,
link-local and null are rejected; a multi-homed host is warned about, naming the alternatives, because
the default route may not be the interface the compute fabric reaches the driver on.

Explicitly **not** address sources: cloud instance metadata (on EC2 the default route already returns
the right address, and off-cloud an unbounded IMDS call hangs) and the Kubernetes downward API (inside
a pod the lookup already returns the pod IP; the API adds a manifest requirement for the same value).
Do not add either back.

R6 is deliberately **not** split per executor: grid, Kubernetes and cloud batch resolve to the same
address by the same reasoning, so separate rungs would be branches with no different behaviour behind
them — and in-cluster or cloud-membership probes are behaviour gated on a signal that cannot be
tested. An earlier draft carried seventeen rungs, each describing a real topology and each arriving
from an adversarial review rather than an observed failure; they were cut.

> **Vestigial code from that cut.** `AgentRpcHostResolver`'s class javadoc still describes a ten-rung
> ladder — a containerized-driver rung, a rootless-podman check, in-cluster Kubernetes detection, a
> cloud-membership probe — and the class still carries the `containerizedValue` / `rootlessValue` /
> `cloudValue` fields, a `DEFAULT_BRIDGE_PREFIX` map and an `AbstractGridExecutor` import for them.
> None of it is reachable: `Probes` exposes only `outboundAddress()` and `interfaceAddresses()`, and
> `ladder()` implements the six rows above. The test file likewise keeps an empty `R3 — the
> containerized driver` section header with no test under it. Read the ladder, not the javadoc.

The advertised endpoint is a **single** `host:port`, not a candidate list. A candidate list would
require withholding the capability token until the fingerprint pin verifies, so that trying N addresses
never offers a bearer token to N hosts. Instead the proxy fails immediately when it cannot connect,
naming the endpoint, and the driver logs the chosen address **with its source label**, so a wrong
answer is one line rather than an hour of silence.

`podman` is unverified: its VM would not boot on the measurement host, so `host.containers.internal`
rests on documentation.

## 7. Tool dispatch and scalability

The stream layer is not the bottleneck. gRPC multiplexes connections on a shared event loop, so an
agent waiting on the model holds an HTTP/2 stream and **no JVM thread**; per-invocation state is small
and every write funnels through a sink that serializes `onNext` and guards post-terminal sends.

**Dispatch is where a thread is consumed.** `ToolDispatcher.call` returns a `String`, so the pool
thread parks for the entire lifetime of the underlying Nextflow task — minutes to hours for a real
module.

| Load | Cost |
|---|---|
| agent waiting on the model | ~0 threads |
| agent inside a tool call | **1 parked thread**, for the task's whole duration |
| 200 such agents | ~200 parked threads, unbounded and growing |

Four defects, ranked: an **unbounded dispatch pool** absorbing load by creating threads until the JVM
degrades, with no push-back and no error naming the cause; the **blocking dispatch SPI**, which is the
ceiling itself; **no global admission control**; and **unbounded diagnostics plus default gRPC
limits**. Practical ceiling: tens of concurrently tool-calling agents, degrading by thread exhaustion
rather than a clean error, against a stated target of hundreds of concurrent invocations and thousands
of idle streams.

The plan is deliberately ordered so the cheap half ships alone.

**M1 — make the ceiling visible and survivable** (no SPI change): a fixed pool plus a bounded queue
sized from `agent.rpc.maxConcurrentTools`; on overflow return a dispatch-level `{"error": …}` tool
result, which the model can already recover from, rather than dropping a frame or growing without
limit; a bounded executor and explicit `maxConcurrentCallsPerConnection` and keepalive on the gRPC
server; rate-limited trace/stderr frames that **report** what they dropped; and a high-water mark
logged at shutdown. M1 does not raise the ceiling — it converts silent thread explosion into
observable queuing.

**M2 — remove the ceiling**: widen the SPI to `CompletableFuture<String> callAsync(…)` with `call` as a
blocking default, so existing runners are unaffected, and have dispatch attach a completion callback
that sends `tool_result` from the completing thread. No thread is parked while a tool task runs, and
concurrency is bounded by the executor's own limits — where it belongs.

Bounding agent concurrency is also a **prerequisite for TLS to scale**: every agent holds an open bidi
stream for its whole life, now each with TLS state and a handshake, so a wide `map` step is what breaks
the broker first, not the crypto.

Non-goals: sharding the broker across processes; a relay transport (deferred for driver reachability,
not for load); parallel tool calls *within* one agent, which is a correlation question; and reconnect
or resume of a broken stream.

## 8. Known limitations

- **A network fault is not recoverable.** The agent holds one stream for its whole life, so an
  interruption fails the task — possibly after paying for model calls — and `errorStrategy 'retry'`
  re-runs the entire loop. Proper reconnect needs idempotent tool replay; the call-ID dedup set is
  groundwork, not a solution. Documented rather than half-built.
- **Completion is stdout-decoded.** A `complete` frame does not by itself mark the task complete today;
  the acknowledged-terminal-state design (§3) is not implemented.
- **A relay transport is deferred**, for drivers that are not reachable from workers at all.
- **`agent.rpc.remoteHost` must still be set** for a driver on a different network from its tasks, a
  Kubernetes driver outside the cluster its pods run in, a cloud-batch driver in a different VPC from
  the compute environment, and a multi-homed submit node whose default route is not the compute fabric.
- **Cacheable local external agents** do not fully represent the packaged runtime and launch command in
  their task identity, so a runtime upgrade can in principle replay an older result. The image tag
  tracking the plugin version narrows this, but does not close it.
- **The image must carry an `agent-rpc` build that understands the current transport flags.** The
  driver always passes `--fingerprint`, or `--insecure` when TLS is off; an older proxy fails with
  `flag provided but not defined`. Turning TLS off is not a workaround.
