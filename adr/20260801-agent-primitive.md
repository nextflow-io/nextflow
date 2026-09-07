# Agent Primitive

- Authors: Paolo Di Tommaso, Ben Sherman
- Status: accepted
- Date: 2026-08-01
- Tags: dsl, agent, llm, language

## Summary

Extend Nextflow's reproducibility model to cover **non-deterministic** compute, by adding a top-level `agent` construct to the DSL: a task-shaped primitive whose body is an agentic tool-calling loop, where each tool call is dispatched as a real Nextflow module invocation inside the same session. The agentic model -- skills, tool invocation, goals, iteration -- then coexists in one artifact with containerization, compute abstraction, portability and caching, rather than living in a harness outside them.

## Problem Statement

Genomics is entering the era of **agentic genomics**: multi-step analyses delegated to autonomous agents that choose their own tools, adapt to intermediate results, and are driven by natural language instead of a hand-written DAG, constrained by libraries of domain-specific *skills* ([Cell Genomics **6**, 101305, 2026](https://doi.org/10.1016/j.xgen.2026.101305)). When the cost of *producing* an analysis collapses, the binding constraint becomes *trusting* it -- the bottleneck moves from pipeline construction to validation.

That makes agentic genomics, first and foremost, a **reproducibility** problem. An agent's tool choices, parameters, iteration count and stopping condition are decided at runtime by a model, so the very things a reproducible pipeline pins down are the things now left open. The surrounding infrastructure that would make such analyses trustworthy -- versioned and discoverable skill registries, deterministic replay of an agent's decisions, signed audit trails, enforced validation tiers, benchmarking, bias detection, equity-aware defaults, viable local-first deployment -- is largely still missing.

**This is the problem Nextflow already solved once, for deterministic compute.** Reproducible, portable analysis was never a property of the science; it was engineered -- containerized tool environments, an executor abstraction that decouples a pipeline from the compute it runs on, explicitly declared inputs and outputs, and content-addressed caching that ties a result to the exact code and data that produced it. Nearly every requirement agentic genomics now raises is that same class of problem, restated for **non-deterministic compute**: replay is caching, audit is lineage, skill distribution is module distribution, local-first deployment is portability.

There is a second convergence, on **scale**. A single agent with one long context does not survive a real workload; the pattern the field is settling on is **agentic map-reduce** -- plan, shard deterministically, map agents in parallel over bounded shards, then reduce -- because it gives coverage by construction and keeps cost proportional to the relevant work rather than to the size of the input ([Devin, *Agentic MapReduce*](https://devin.ai/blog/agentic-map-reduce)). Its guiding rule is one a workflow engineer will recognize immediately:

> Put agents where reasoning is required -- synthesizing the decomposition function, inspecting the shards, and the reduction. Everything else is deterministic.

That is fan-out/fan-in with reasoning at selected nodes, which is a description of a dataflow graph. And the fit is deeper than the shape: in a dataflow model a node *awaits its inputs and fires when values arrive*, which is exactly what an agent is -- something that waits for a signal to process, runs, and emits. Nextflow already provides that execution model, along with the parallelism, backpressure, fan-in and ordering semantics that come with it, so an agentic map-reduce needs no orchestration layer of its own: the shard is a process, the map is an agent over a channel, the reduce is the same agent over a collected one. Expressing this pattern outside the engine means rebuilding concurrency control, work distribution and gather semantics that a dataflow runtime already has.

A third convergence, on **packaging**. An agent's durable value is not a prompt someone typed once; it is a reviewed unit -- an instruction, a model, an output contract, the skills that constrain it and the tools it may call -- and that unit is only reusable if those parts travel together, versioned, so a second lab can run it and get the same behaviour. The skill-registry requirement above is, at bottom, a packaging and distribution problem: units that are findable, versioned, tested, documented and carrying enough metadata to be selected by something other than their author.

This is a solved problem in Nextflow, and solved in the same shape the paper is asking for. Modules already encapsulate a unit of analysis with its dependencies and metadata; `include { x as y }` composes and aliases them; the registry and semantic versioning resolve them; config selectors parameterize them from outside without editing them; and nf-core has demonstrated the community governance the paper explicitly proposes as the model for skill registries. Extending that to agentic compute means an agent is a **module** too -- encapsulated together with its skills and its tools, included under an alias, versioned and resolved like everything else -- rather than a per-project prompt file with a bespoke loader.

The overlap is therefore not incidental, and the gap is one of placement rather than capability. Today agentic compute sits *outside* the engine: an external harness calls in over MCP or a platform API to launch pipelines, holding its own provenance model, its own retry semantics, its own concurrency model, and no shared work directory -- so every requirement above has to be solved again inside that harness, where none of the existing machinery reaches.

**This ADR proposes the inverse: make the agent a node in the workflow graph.** Rather than build a parallel stack for agentic compute, extend the reproducibility framework to cover non-deterministic steps, so that skills, tool invocation, goals and iteration coexist in a single artifact with containerization, compute abstraction, portability and caching -- and an agent's decisions become subject to the same execution, caching, provenance and container machinery as every other task.

## The infrastructure agenda, mapped

The claim above -- that the overlap is real and the gap is placement -- is only worth making if it survives contact with the specifics. Below are the infrastructure requirements the perspective identifies, against what Nextflow already provides and what this design does with it. "Substrate" means the mechanism exists and is load-bearing today; "gap" means this ADR does not address it and something else must.

| Requirement | Nextflow substrate | This design |
|---|---|---|
| **Skill registries** -- standardized skill metadata, semantic discoverability, versioning; GA4GH and FAIR alignment; DOI-minted releases. The paper explicitly proposes modelling governance on nf-core: *"agent-native skill registries should adopt analogous governance structures"* | nf-core's curation model, and the module registry with versioned modules, `meta.yml` metadata and semantic resolution ([module system ADR](20251114-module-system.md)) | A tool **is** a registry module -- no parallel packaging format. Tool descriptions and I/O schemas come from `meta.yml` / registry `ModuleMetadata`. Skills are versioned module-local `SKILL.md` bundles. DOI minting and formal GA4GH metadata are registry work, not language work |
| **Deterministic replay** -- *"given a logged sequence of decisions, it must be possible to reproduce the exact execution path, even if the original selection was mediated by an agent"*; required for regulatory audit, debugging and multi-site concordance | `-resume` over a content-addressed task hash; the work directory as the durable record of every step | **Delivered, and the sharpest fit in the paper.** Because an agent invocation is a task, its canonical identity -- runner, model, instruction, goal, prompt, output schema, skill-content fingerprint, tool fingerprint -- is folded into the hash, and a replay is served only for that exact configuration. Editing a tool invalidates it |
| **Signed reproducibility bundles and audit trails** -- BioCompute Objects, version-locked dependencies | Lineage records ([data lineage ADR](20250508-data-lineage.md)); container digests; the work dir | Substrate: an agent run emits an `AgentRun` lineage record naming the resolved model, tools and skills, and tool tasks are container-pinned. BCO emission and bundle signing are **gaps** |
| **Enforced tier boundaries** -- research / benchmarked / clinical grade, *"preventing agents from applying research-grade skills in clinical contexts without explicit override"* | Config scopes and selectors; module metadata as the place a tier would be declared | **Gap.** Nothing declares or enforces a tier. This is the most consequential unaddressed item, and module metadata is where it belongs |
| **Formal benchmarking infrastructure** -- standardized datasets, confidence intervals, controlled multi-site evaluation; plus held-out partitions and rotated benchmarks to avoid circular tuning | `nf-test`; nf-core's community benchmarking practice; portability across sites, which is what makes multi-site evaluation mechanically possible at all | **Gap** as infrastructure. The primitive contributes only the precondition: because a tool is a tested module and an agent result is an ordinary channel value, an agent is benchmarkable by the same harness as a pipeline |
| **Out-of-distribution and bias detectors as first-class components** -- population mismatch, sequencing-platform drift, tissue-context mismatch; the conceptual machinery (calibration curves, conformal prediction, density-ratio estimators) exists and *"what is missing is their routine integration into agentic pipelines"* | Processes and channels: a detector is just another module in the graph | Substrate. Because agent output is a typed channel value, a detector composes downstream of an agent exactly as it would after a caller. What is missing is the **library** of such detectors and any obligation to run them -- a community/registry concern |
| **Equity as a systems requirement** -- equity-aware defaults, population-aware model selection with population metadata, bias monitoring during execution, measurable metrics such as the proportion of skills validated on non-European populations | Module metadata; declared, reviewable tool sets | Partial, and indirect: because the tool set is *declared* rather than discovered, an agent cannot silently reach for the most abundant reference resource -- it selects from a curated set a human reviewed. Population metadata and equity metrics themselves are **gaps** |
| **Reproducibility by design, and lightweight preregistration** -- a plan expressed as a fixed sequence of skill invocations deposited at a registry (e.g. OSF), with deviations *"automatically detectable as differences between the registered and executed skill traces"* | Pipelines are already versioned, executable, git-resolved artifacts; the DAG and lineage are the trace | Strong fit, unclaimed. A Nextflow script with pinned module versions *is* a preregisterable plan, and lineage is the executed trace, so the registered-vs-executed diff the paper wants is mechanically available. Nothing implements the comparison |
| **Model hosting for low-resource settings** -- intermittent connectivity, documented minimum hardware, pre-packaged equity-aware defaults, local-first computation so *"genetic data never leave the researcher's machine"* | The executor abstraction and containers: the same pipeline runs on a laptop and on a cluster, which is precisely the local-first property being asked for | Substrate, plus one addition: the `AgentRunner` SPI makes the model endpoint swappable like an executor, so the same pipeline can target a hosted API or a local model without touching pipeline code |

Two conclusions follow. First, Nextflow is not incidental to this agenda -- for replay, audit, portability, local-first execution and preregistration, the mechanism the paper asks for already exists and is in production. Second, the items that remain gaps cluster in one place: **skill and module metadata** (tiers, population scope, equity metrics) and **community governance** (benchmark datasets, detector libraries). Those are registry and ecosystem problems. They are not reasons to keep the agent outside the engine; they are the next ADRs.

## Why a primitive and not a process?

A Nextflow pipeline could implement agentic workflows today, using processes that call an agent harness via CLI or API. However, this approach requires a significant amount of boilerplate to implement while preserving Nextflow best practices:

- Supporting one or more agent harnesses in a portable and scalable manner
- Designing the agent-in-a-process to cache only the relevant inputs
- Enforcing the declared output structure with a JSON schema
- Recording an agent run correctly as a `TaskRun` lineage record
- Enabling dynamic workflows over a set of modules

This "boilerplate" is significant enough to warrant a language primitive, as implementing all of this pipeline code would distract from the pipeline's core purpose of defining the data flow.

An agent primitive solves all of these issues in the runtime:

- An `agent` definition is more concise than an equivalent `process` definition. It provides built-in caching and output validation.
- The agent runner (harness) is a plugin extension, allowing it to be configured independently from the pipeline code (like executors).
- Agent runs are recorded in lineage with a dedicated `AgentRun` lineage record.
- Agents can run modules in the same compute environment as the main run.

The following examples highlight agentic use cases that would be difficult to implement without an agent primitive:

- `examples/agents/11_contig-filter`
- `examples/agents/12_isolate-triage`
- `examples/agents/15_map-reduce`

## Goals or Decision Drivers

- **Extend the reproducibility framework, do not fork it**: a non-deterministic step must inherit the same guarantees a deterministic one gets -- container-pinned dependencies, an executor-agnostic work dir, a content-addressed cache key, and a lineage record.
- **First-class language construct**: agents compose with processes and other agents through the existing channel/workflow model.
- **Task-shaped primitive**: one invocation per input record, executed as a `TaskRun`, so parallelism, caching, resume, retry, lineage and reporting are inherited rather than rebuilt.
- **Agentic map-reduce with no orchestration layer**: an agent is a dataflow node awaiting its inputs, so fan-out over a channel, parallel map and `collect()`-style fan-in come from the execution model rather than from a scheduler written on top of it.
- **Agents and skills are modules**: an agent is distributable as a versioned module carrying its own skills and tools, includable and aliasable like a process, and configurable from outside via the `agent` scope and selectors -- reusing module resolution rather than inventing a skill format.
- **Tool calls are real dataflow nodes**: the module runs through the standard executor/container/retry/cache machinery and its output flows back to the agent.
- **Modules are tools**: a module's description, inputs, and outputs come from its `meta.yml` / registry `ModuleMetadata`, so there is no parallel tool-metadata format to maintain.
- **Runner behind an SPI**: the agent harness is swappable like an executor; core does not depend on a specific harness.
- **Backward compatibility**: zero impact on existing `process`/`workflow`/`function` declarations.

## Non-goals (v1)

- Long-lived conversational agents and cross-invocation memory (one invocation → one result).
- Channel-aware orchestrator agents (an agent does not subscribe to channels mid-loop).
- Agent-to-agent invocation as tools (composition is via channels).
- Cost tracking, token budgets, prompt analytics.
- Streaming partial outputs.
- The registry- and community-side items in the agenda above: validation-tier declaration and enforcement, benchmarking infrastructure, OOD/bias detector libraries, population and equity metadata, BCO bundle emission, preregistration diffing. These are metadata and governance work that this language change should enable, not absorb.

## Solution

A new `agent` keyword, lowered to a V2 `TaskProcessor` whose body is the agent loop; declared `tools` are namespaced references -- `nf:module_run[:<NAME>]` over the processes in scope, whether declared locally or `include`d from a local or registry module, plus `fs:`/`shell:` for the runner's own tools -- and the selected modules are pre-wired as dataflow gateways; the harness lives behind the `AgentRunner` SPI in a plugin.

```nextflow
agent shouty {
    model 'openai/gpt-5-mini'
    instruction 'Call the `uppercase` tool, then reply with only its result.'
    tools 'nf:module_run'

    input:
    request: String

    output:
    answer: String

    prompt:
    "${request}"
}
```

## Rationale & discussion

The defining choice is **agent-as-task**. Each invocation is one record-in / result-out node; multi-turn reasoning happens inside it, and multi-agent pipelines compose like multi-process pipelines. This reuses channel composition and adds exactly one node type -- and, crucially, it is what makes the reproducibility machinery apply to agent decisions at all.

Modules-as-tools over a bespoke tool format reuses the `meta.yml` schema already validated by `ModuleSchemaValidator`, so tool descriptors come for free. This is the same interoperability argument the paper makes from the other direction -- *"a skill can wrap a Nextflow module […] the agent paradigm must build on, not displace, two decades of community investment in reproducible infrastructure"* -- except that here the module **is** the skill, with no wrapper layer.

### Does the primitive satisfy the definition?

The infrastructure agenda is mapped above; this is the narrower question of whether the construct is the right *shape*. The paper's definition is four necessary conditions, and it makes them empirically testable via a **perturbation test**: a system fails to qualify as agentic if, given identical input with a perturbed intermediate result, it does not alter its execution strategy.

| Condition | Mechanism here |
|---|---|
| **Autonomy** -- decisions during execution, not a static specification | The body is a tool-calling loop: which tool runs, with which arguments, and how many turns, are decided at runtime up to `maxIterations` |
| **Domain constraint** -- a structured library of validated operations, not code generated ad hoc | `tools` / `skills` resolve to in-scope processes and versioned registry modules; the agent cannot invent a tool outside the declared set. Enforced by the language, not by prompt discipline |
| **Iterative refinement** -- evaluate intermediate results, diagnose, self-repair | Tool results, including dispatch-level errors, are returned to the model so it can adapt within the loop -- which is what makes the perturbation test pass rather than fail |
| **Natural-language mediation** | `prompt:` with typed interpolation, plus `instruction` / `goal` |

So the construct qualifies on the paper's own terms, while remaining a task -- which is the whole point.

It is worth stating the counter-position plainly, because the perspective makes it directly:

> Nextflow, Snakemake and Galaxy provide powerful, reproducible pipeline orchestration, but their workflows are specified in advance by a human developer. […] Agentic genomics builds on workflow infrastructure (an agent may invoke a Nextflow pipeline as one step in a larger analysis) but is not reducible to it.

That is correct about workflow managers *as they are*, and it is the reason this is a language change rather than a library: what makes the DAG insufficient is precisely that it is fixed in advance, so the fix is a node type whose edges are chosen at runtime. Notably, the paper treats deterministic replay as an unmet gap across all four systems it surveys (CellAtria, AutoBA, Bio-Copilot, ClawBio), every one of which sits *outside* a workflow engine. The claim here is that this is not a coincidence: putting the agent *inside* the engine is what makes replay fall out of existing machinery instead of having to be built.

## Links

- Agentic genomics framing: Corpas, M., Guio, H., and Fatumo, S. (2026). *Agentic genomics: From pipeline automation to autonomous validation.* Cell Genomics **6**, 101305. <https://doi.org/10.1016/j.xgen.2026.101305>
- Agentic map-reduce as the scaling pattern: [Devin, *Agentic MapReduce*](https://devin.ai/blog/agentic-map-reduce)
- User guide: [`docs/agent.mdx`](../docs/agent.mdx)
- Design: [`adr/specs/agent-design.md`](specs/agent-design.md) -- the primitive: lowering to a task, typed I/O, tools, skills, modules, configuration, caching and lineage
- Design: [`adr/specs/agent-runners.md`](specs/agent-runners.md) -- the `AgentRunner` SPI, runner selection, and the two shipped harnesses
- Design: [`adr/specs/agent-rpc.md`](specs/agent-rpc.md) -- canonical task execution, the driver broker protocol, transport security and scalability
- Examples: [`examples/agents/`](../examples/agents) (`validate.sh` runs them end to end; `-r` also checks that each replays from cache on `-resume`)
- Related ADRs: [`record types`](20260306-record-types.md), [`module system`](20251114-module-system.md), [`data lineage`](20250508-data-lineage.md), [`type system`](20260501-type-system.md)
- Inverse pattern this design inverts: [`colbyford/nf-foundry-workflow`](https://github.com/colbyford/nf-foundry-workflow) -- Foundry agents calling Seqera MCP to launch Nextflow pipelines
- Runtime engine: [`langchain4j`](https://docs.langchain4j.dev/tutorials/agents)
