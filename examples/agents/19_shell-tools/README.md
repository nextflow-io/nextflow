# 19_shell-tools — runner-native file and shell tools

```groovy
tools 'fs:*', 'shell:bash'
```

## What it shows

Every other tool example brokers its tools back to the driver: `nf:module_run` marshals the model's
arguments into channel values and runs a real Nextflow task. These tools do the opposite — they are
the **runner's own**, activated by name and executed where the agent already is. On `pi` that is the
SDK's builtins inside the agent container, so no tool call crosses the RPC link.

The task is arithmetic no model can do by reading. `data/contigs.fa` holds 8 contigs and 2,670
bases; GC content and N50 are exact quantities over every base. Without a shell the agent can only
estimate, and an estimate here is just a wrong answer. With `awk` it computes.

That is the case the tools redesign was argued from: *"grep / sed / awk a file to find targeted
information, or write custom code to compute things it can't do directly with tokens."*

## How the FASTA reaches the agent

The agent declares `contigs: Path`, which means for an agent exactly what it means for a process:
the file is **staged into the agent's task directory** and bind-mounted into the runner container,
and both the prompt interpolation and the input JSON render it as the plain name `contigs.fa`. So
the agent's shell opens it directly — no `write` step, no copy of the sequence through the context
window.

Two consequences worth carrying:

- Beyond its declared inputs, `fs:*` and `shell:bash` on `pi` see **only what the agent itself
  creates** in its sandbox. They are not a general window onto pipeline data.
- Module outputs reach the model differently again: small text/JSON results are **inlined into the
  tool result** (`ToolOutputReader`, 32 kB cap, `.json`/`.tsv`/`.txt`/…), so an agent gets those
  numbers without ever opening a file. Anything else arrives as an opaque path handle.

On `langchain4j` the picture differs — the loop runs in the driver JVM, so `fs:*` reaches the real
filesystem behind `SandboxGuard` (the work dir, the agent's staged inputs, and module outputs).
`shell:bash` is not available there at all.

## Expected result

The numbers are fixed — `data/contigs.fa` is generated from a pinned seed — so this example is
checkable rather than merely plausible:

```
contigs   = 8
bases     = 2670
GC%       = 55.51
N50       = 620
longest   = 800
```

Verify independently:

```bash
grep -c '^>' data/contigs.fa
awk '/^>/{next}{n=length($0); tot+=n; gc+=gsub(/[GCgc]/,"")}END{printf "%d bases, %.2f%% GC\n", tot, gc*100/tot}' data/contigs.fa
```

A wrong N50 is the interesting failure: it is the one statistic needing a sort and a running total
rather than a single pass, so it is where an agent that quietly stopped using the shell gives itself
away.

## Why `pi` only

`shell:bash` is rejected at agent-build time on `langchain4j`. That loop runs in the **driver JVM**,
so a shell tool there would execute model-authored commands on the driver host with no container
boundary — and on Kubernetes or Batch the driver pod does not have the pipeline's tooling anyway. In
the `pi` container, the container *is* the sandbox.

| runner | `fs:*` bounded by | `shell:bash` |
|---|---|---|
| `pi` | the session `cwd` plus the container | the container only |
| `langchain4j` | `SandboxGuard` (work dir + staged inputs + module outputs) | not available |

## Run

```bash
export OPENAI_API_KEY=...
nextflow run .
```

Requires the `pi` runner image the `nf-agent-pi` plugin declares and a running
container engine.
