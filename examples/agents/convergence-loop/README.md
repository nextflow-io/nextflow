# convergence-loop — a true convergence loop (same tool, many runs)

An agent that re-runs the **same** tool many times, varying one parameter and
reading the metric each time, until it converges on the optimum. Fully offline.

## Purpose / What it demonstrates

This is a *loop* in the strict sense — not tool composition (different tools
chained once, as in [`goal-directed`](../goal-directed)) but **iterated execution of one
module to optimise a result**. The agent runs the `score_threshold` tool over
and over, each time with a different `threshold`, reads the F1 it returns, and
steers toward the value that maximises F1: a genuine coarse-scan-then-refine
search where each next guess depends on the previous result.

The key enabling fact: **iterative tuning needs no core change today** because
the tunable knob is a *declared process input*. `module_run` feeds declared
inputs per call, so the agent can vary `threshold` on every invocation. (Tuning
a *stock* nf-core module's hidden `task.ext.args` would instead require a planned
core feature — surfacing `ext.args` in the tool schema — and is out of scope
here.)

The task is a real one: choosing the score cutoff that maximises the **F1** of a
binary filter against a truth-labelled set — a standard "operating point"
decision (e.g. a variant or quality-score threshold). The F1-vs-threshold curve
has a single interior peak the model cannot know without running the tool, so it
*must* iterate to find it.

## How it works

1. **The `score_threshold` process** reads a fixed, committed labelled dataset
   (`scores.txt`, columns `score label`) and computes precision / recall / F1 of
   the rule `score >= threshold`, returning a one-line string. It is a local
   `exec:` process — no container.

2. **Fixed dataset, single knob.** The dataset path is a config value
   (`params.scores`), **not** an agent input, so every tool call is just
   `{"threshold": <x>}`. That keeps the loop a clean single-parameter search and
   removes any chance of the model mangling a long constant path across dozens of
   calls. (A local process tool supports only scalar inputs, so the knob is a
   `BigDecimal` — the JSON number type.)

3. **The `tuner` agent** has a role-only `instruction` ("the tool is the only
   way to learn F1 — never guess; scan coarsely then refine") and a `goal` ("find
   the threshold that maximises F1; converge on the best; report it"), with
   `maxIterations 30`.

4. **Verified convergence:** in a real run the model made **32 calls** — a
   coarse scan `0.0 → 1.0`, then refinement around `0.45–0.56`, then fine steps
   `0.472 … 0.484` — converging on **threshold 0.48, F1 = 0.8889** (the true
   optimum), with zero failures.

## Key concepts

| Concept | In this example |
|---|---|
| Convergence loop | The same tool run many times, converging on an optimum |
| Loop ≠ composition | Iteration of one tool vs chaining different tools once |
| Declared-input knob | `threshold` is a process input ⇒ tunable per call, no core change |
| Fixed dataset | `params.scores` (not an agent input) ⇒ tool call is just `{threshold}` |
| Scalar-only local tool | Knob is `BigDecimal`; result is a `String` |
| The `ext.args` limit | Stock-module hidden params would need a core feature (see `docs/agent.mdx`) |

## Running it

**Requirements:** an OpenAI API key and the `nf-agent` plugin. No container
runtime and no data fetch — `scores.txt` is committed alongside the example.

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected output (the model's report):

```
OPTIMAL=Optimal threshold = 0.48 — precision 0.80, recall 1.00, F1 0.8889
```

See [examples/agents/README.md](../README.md) for the dev-build (run-from-repo)
instructions.
