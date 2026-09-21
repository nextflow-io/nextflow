# fan-in-parity — an `agent` reduces exactly like a `process`

Side-by-side proof that a tool-free `agent` inherits **canonical Nextflow cardinality**.
A `process` reducer and an `agent` reducer consume the **same** collected channel — a
`Bag<Finding>` produced by `collect()` — and **both fire exactly once** over the whole
bag. Nothing about the fan-in is agent-specific; only the body differs (deterministic
Groovy vs an LLM call).

## The two reducers (identical input contract)

```nextflow
process reduce_with_process {          agent reduce_with_agent {
    input:                                 model 'openai/gpt-5-mini'
        findings: Bag<Finding>             instruction '...'
    output:                                input:
        report: String                         findings: Bag<Finding>
    exec:                                  output:
        report = "...${findings.size()}"       report: String
}                                          prompt: "...${findings.size()}..."
                                       }
```

```nextflow
def bag = channel.of(f1, f2, f3).collect()   // one value item: a Bag of all 3
reduce_with_process(bag)                       // fires once
reduce_with_agent(bag)                         // fires once
```

## Why they behave the same

This is the ordinary process rule, not a special "reduce" primitive:

- `collect()` turns a queue channel into a **value (singleton) channel** holding one
  `Bag` of all items.
- A process fires N times = the length of its longest **queue** input; **value** inputs
  are broadcast (reused), not consumed. All-value inputs ⇒ **one** invocation.
- The agent reuses the exact same input-matching (`ProcessInputsDef.isSingleton()`), so a
  `collect()`/value input fires it once — identical to the process.

The `Bag<E>` type is not agent-specific either: `Channel.collect()` is typed
`Value<Bag<E>>`, so a typed *process* reducing over `collect()` declares `Bag<E>` too.
(Order over the bag is not guaranteed and is hashed order-independently — again, standard
Nextflow, so resume survives reordering.)

## Observed (verified end-to-end)

Fresh run — each reducer is one task:

```
[c4/4bf198] Submitted process > reduce_with_process
[aa/805a26] Submitted process > reduce_with_agent
PROCESS reducer => combined 3 findings [f1, f2, f3]
AGENT   reducer => The analysis shows low contamination, minimal adapter content, and good coverage.
[SUCCESS] completed=2 failed=0 cached=0
```

`-resume` — both cache identically:

```
[c4/4bf198] Cached process > reduce_with_process
[aa/805a26] Cached process > reduce_with_agent
(Submitted=0  Cached=2)
```

## Running it

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf            # completed=2 (one process task + one agent task)
nextflow run main.nf -resume    # cached=2
```

Requires the `nf-agent-pi` plugin (in `nextflow.config`) and an OpenAI key for the agent
reducer; the process reducer is a plain `exec:` and needs neither.
