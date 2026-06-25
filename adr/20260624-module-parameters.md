# Module Parameters and Tool Arguments

- Authors: Paolo Di Tommaso
- Status: proposed — pending decision between Option A and Option B (author recommends B)
- Date: 2026-06-24
- Tags: modules, parameters, tools, arguments, configuration, records
- Related: [Module System](20251114-module-system.md), [Workflow params](20250825-workflow-params.md), [Record types](20260306-record-types.md)

## Summary

Unify two previously separate drafts — *Module Parameters* ([PR #6769](https://github.com/nextflow-io/nextflow/pull/6769)) and *Tools Arguments* (an unpublished draft, not committed to `adr/`) — into a single mechanism: a tool flag is just a parameter whose string form happens to be `-K <value>`. Both drafts agree on the host (`params {}` block), the override binding, and the custom-type machinery; they differ only in **how grouping/nesting is expressed**. This ADR presents the two candidate designs as **Option A** (a `group` attribute + a single `ToolArg` type) and **Option B** (record-typed params with `ToolOpt`/`ToolArgs`), then compares them. This ADR supersedes both drafts.

## Problem Statement

Two adjacent problems were being solved by two separate designs:

1. **Module parameters** — modules need a documented, type-safe way to declare and override configurable values, replacing scattered `params.x` assignments.

2. **Tool arguments** — modules wrap one or more command-line tools (e.g. `bwa`, `samtools`), and their CLI flags are today passed through `task.ext.args`/`ext.args2` as raw strings. This mechanism has fundamental shortcomings:
   - **Buried and scattered.** The values live in `ext.args` assignments spread across several config files (`nextflow.config`, profiles, institutional configs, `conf/modules.config`), with no single place that declares what a module accepts.
   - **Undocumented.** A raw string like `-K 100000000 -Y` carries no description of what each flag does, its type, or its allowed values.
   - **No schema.** There is no machine-readable definition of a module's arguments, so they cannot be validated, introspected, or surfaced by tooling and IDEs.
   - **Opaque semantics.** Arguments are concatenated strings — easy to mistype, impossible to type-check, and meaningless to anything but the underlying tool.

Keeping module params and tool args as two parallel concepts (`params` vs `tools.<tool>.args`) compounds this: two specification schemas, two override syntaxes, and two mental models for the same thing — *a named, typed, overridable input to a module*.

**Goal:** expose tool arguments (and module settings generally) as **first-class parameters** — declared in one place, documented and backed by a schema for validation, and accessible programmatically and consistently across the CLI, config, and APIs.

## Decision Drivers

- **Replace `task.ext.args`.** Superseding the opaque `ext.args`/`ext.args2` strings with a documented, typed, introspectable mechanism is the primary goal — not merely offering an alternative.
- **One concept, one syntax.** A tool argument is a parameter with a particular serialization. Don't invent a second namespace for it.
- **Reuse the `params {}` block.** The typed `params {}` block already exists (workflow params ADR). Module params should be the same block.
- **Namespacing.** Modules wrap multiple tools (`bwa` + `samtools`) and have logically related params; both need grouping to avoid name collisions and document tool boundaries.
- **Minimize new grammar / concepts.** Prefer mechanisms the language already has.
- **Unambiguous overrides.** Config/CLI overrides must address the right module unambiguously.
- **Extensibility.** Tool-argument types should be instances of an open custom-type mechanism, not hard-coded specials.

## Non-goals

- **Immediately removing `task.ext.args`.** Replacing it is a goal (see Decision Drivers), but existing pipelines keep working through a deprecation period — this ADR does not delete it outright.
- **Changing the workflow-level `params {}` block** ([Workflow params](20250825-workflow-params.md)) — this ADR extends its use to modules, it does not alter workflow params.
- **Replacing `nextflow_schema.json`.**
- **Resolving the subcommand-collision, repeated-flag, and custom-type-API questions** — these are tracked under [Open Questions](#open-questions).

---

# Common foundations (shared by both options)

These elements are identical regardless of which grouping design is chosen.

### Host: the `params {}` block

Module parameters are declared in a module-level `params {}` block, reusing the workflow-params syntax (`name: Type` / `name: Type = default`). The two options below differ only in how grouped/tool params are declared *inside* this block.

> The original *Module Parameters* draft (PR #6769) weighed a **process-level** `params:` block against this **module-level** `params {}` block, and the module-level block is adopted here for continuity with workflow params. Note that the "Option A / Option B" labels in *this* ADR denote a **different axis** — how grouping is expressed — than the A/B options in PR #6769 (which were about host placement).

### Specification attributes (`meta.yaml`)

| Attribute | Required | Applies to | Description |
|-----------|----------|------------|-------------|
| `name` | Yes | all | Parameter / field / tool-option name |
| `type` | No (default `string`) | all | `boolean`, `integer`, `float`, `string`, `file`, `path`; `record` or a record type; `ToolArg` (A) / `ToolOpt`,`ToolArgs` (B); or a custom type |
| `description` | No | all | Human-readable docs |
| `example` | No | all | Example value (docs / IDE) |
| `default` | No | all | Value used when unset |
| `enum` | No | scalar / tool option | Allowed values (validated) |
| `group` | No | all (**Option A only**) | Namespace → `params.<group>.<name>` |
| `fields` | When record-typed | record / `ToolArgs` (**Option B only**) | Nested member declarations |
| `<Type>:<key>` | No | custom types | Type-specific attribute, e.g. `ToolArg:prefix` / `ToolOpt:prefix` |

`description`, `example`, and `enum` are documentation/validation metadata that live only in `meta.yaml`; the script `params {}` block carries only `name: Type = default`. The `meta.yaml` spec and the script declarations are expected to agree — how divergence is handled is an [open question](#open-questions).

### Binding resolution

The override **path** (`params.<group>.<name>`) is identical in both options; only the prefix differs by context.

**Within a workflow** the path is scoped by the **name under which the process is included** (resolving the naming ambiguity of a module-level block); for a **standalone `module run`** there is a single implicit process, addressed directly. These overrides use Nextflow's single-dash *runtime options* (`-params.…`, `-process.…`) — **not** the double-dash `--<name>` pipeline-parameter syntax, which would create a flat param literally named `params.bwa.K`.

**Standalone module run:**

```bash
nextflow module run <module/bwa_mem> \
    -params.foo.batch_size=1000 \
    -params.bwa.K=100000000 \
    -params.bwa.Y
```

**Within a workflow** — scoped by the included (possibly aliased) process name:

```groovy
include { BWA_MEM } from './modules/bwa'
include { BWA_MEM as BWA_MEM_2 } from './modules/bwa'
```

```groovy
// config
process {
    withName: 'BWA_MEM'   { params.bwa.K = 100000000; params.bwa.Y = true }
    withName: 'BWA_MEM_2' { params.bwa.K = 50000000 }
}
```

```bash
# CLI
nextflow run main.nf \
    -process.BWA_MEM.params.bwa.K=100000000 \
    -process.BWA_MEM.params.bwa.Y \
    -process.BWA_MEM_2.params.bwa.K=50000000
```

The same module included twice can be configured independently.

### Custom parameter types

Tool-argument types are instances of an open mechanism. A custom parameter type defines two operations:

- **parse**: `String → value` — how a CLI/config string is coerced into the typed value.
- **serialize**: `value → String` — how the value renders when interpolated into a script.

A custom type may also extend the `meta.yaml` param spec with its own attributes, namespaced as `<Type>:<key>` (valid YAML, since a colon not followed by a space is part of the key). Nextflow passes the collected `<Type>:*` keys to that type's logic; unknown namespaces are an error. Both options use this to override a tool option's CLI prefix (`ToolArg:prefix` / `ToolOpt:prefix`).

### Parameter semantics

- **Typing & defaults:** parameters are typed; defaults must be literals/constants or references to other params. CLI strings are coerced per the declared type (or the custom type's `parse`).
- **Required vs optional:** following the [Workflow params](20250825-workflow-params.md) rules, a param with no default is *required* (the run fails if unset); a `Type?` param is *optional* (`null` if unset); a `Boolean` defaults to `false`. **Tool options are optional by nature** — an unset `ToolArg`/`ToolOpt` contributes nothing to the command line and is not treated as a missing required param.
- **Validation & errors:** supplying a param not declared in the `params {}` block or config is an error; values are validated against the declared type and any `enum`, with failures reported at resolution time, before execution (per the workflow-params resolution ordering).
- **Parameter vs input:**

  | Aspect | Parameter | Process input |
  |--------|-----------|---------------|
  | Resolved | Once, before execution | Per task |
  | Mutability | Immutable | Varies per task |
  | Source | Config, CLI, params file, defaults | Channels / values |

  A record-typed *param* (e.g. `ToolArgs`) is configuration; a record-typed process *input* ([Record types](20260306-record-types.md)) is per-task data. They share the record machinery but are distinct roles.

- **Precedence** (lowest → highest): default in `params {}` (or `meta.yaml`) → config → params file → command line.
- **Params file:** nested groups map to nested objects, e.g.
  ```json
  { "bwa": { "K": 100000000, "Y": true }, "output_format": "cram" }
  ```

### Tool option serialization (common to both `ToolArg` and `ToolOpt`)

A single tool-option value, whatever the type is named, renders to its command-line form. The **name is the tool option name**, and the group/record name is the **tool name** (`bwa`, `samtools`) — which is why `params.bwa` renders that tool's full argument string.

| Value | Renders as |
|-------|------------|
| `K = 100000000` | `-K 100000000` |
| `Y = true` | `-Y` |
| `output_fmt = "cram"` | `--output_fmt cram` |
| unset / null / false | `` (omitted) |

**Prefix inference** (overridable via `<Type>:prefix`): a single-character name uses `-`; a name of two or more characters uses `--`. **Boolean** `true` emits only `prefix+name`. An `enum` constrains the allowed values; a `default` supplies the value when unset. Components are accessible via `.name` and `.value`.

### Programmatic access

Within the script, params are read like any other value. Assuming the configuration `params.foo.batch_size = 1000`, `params.bwa.K = 100000000`, `params.bwa.Y = true`, and `output_format` left at its declared default `'bam'`:

| Expression | Kind | Resolves to |
|------------|------|-------------|
| `params.output_format` | scalar param | `"bam"` (declared default) |
| `params.foo.batch_size` | grouped scalar | `1000` (an `Integer`, not a flag) |
| `params.bwa` | tool group | `"-K 100000000 -Y"` (space-joined set options) |
| `params.bwa.K` | tool option | `"-K 100000000"` (string form) |
| `params.bwa.K.name` | option name | `"K"` |
| `params.bwa.K.value` | option value | `100000000` (typed underlying value) |
| `params.bwa.Y` | tool option (boolean) | `"-Y"` |
| `params.bwa.Y.value` | option value | `true` |

So a tool option behaves as its rendered flag in string context, while `.name`/`.value` expose the parts for programmatic logic; a generic grouped param (`params.foo.batch_size`) is just its typed value — only tool options/groups render as flags.

---

# Option A — `group` attribute + `ToolArg`

Grouping is an explicit attribute. Generic params and tool options both declare a `group`; tool options additionally use the built-in `ToolArg` type. (`group` is the *only* namespacing mechanism — there is no `ToolArg:tool` key.)

### meta.yaml

```yaml
params:
  - name: output_format          # ungrouped scalar -> params.output_format
    type: string
    enum: ["sam", "bam", "cram"]
    default: "bam"
    description: "Output format"

  - name: batch_size
    type: integer
    group: 'foo'
    description: "The batch_size description"
    example: 1000
  - name: other_size
    type: integer
    group: 'foo'
    description: "The other_size description"
    example: 2000

  - name: K
    type: ToolArg
    group: 'bwa'
    description: "The bwa -K CLI option"
  - name: Y
    type: ToolArg
    group: 'bwa'
    description: "The bwa -Y CLI option"
  - name: threads
    type: ToolArg
    group: 'samtools'
    ToolArg:prefix: '-@'
    description: "Number of sort threads"
```

`group` produces the nesting level: `params.foo.batch_size`, `params.bwa.K`, `params.samtools.threads`.

### Script definition

A `group` is expressed as a nested block inside `params {}`:

```groovy
params {
    output_format: String = 'bam'      // ungrouped -> params.output_format

    foo {                               // -> params.foo.*
        batch_size: Integer = 1000
        other_size: Integer = 2000
    }

    bwa {                               // -> params.bwa.*
        K: ToolArg
        Y: ToolArg
    }

    samtools {
        threads: ToolArg
    }
}
```

### The `ToolArg` type

`ToolArg` is the single built-in type for a tool option, serializing per the common rules above. A group whose members are `ToolArg` renders, when interpolated whole, as the space-joined serialization of its set members:

```groovy
// params.bwa  ->  "-K 100000000 -Y"
bwa mem ${params.bwa} -t $task.cpus $index $reads
```

(This group-concatenation is an implicit rule for groups of `ToolArg`; individual options remain accessible as `params.bwa.K`.)

---

# Option B — record-typed params (`ToolOpt` + `ToolArgs`)

There is no `group` attribute and no nested-block syntax. The `params {}` block stays **flat**; grouping is achieved by giving a parameter a **record type**, reusing the [Record types](20260306-record-types.md) machinery. Tool options use two built-in types: `ToolOpt` (a single option) and `ToolArgs` (a record of `ToolOpt`).

### meta.yaml

A record-typed param declares its members under `fields:` (same shape as the top-level `params:` list):

```yaml
params:
  - name: output_format          # ungrouped scalar -> params.output_format
    type: string
    enum: ["sam", "bam", "cram"]
    default: "bam"
    description: "Output format"

  - name: foo
    type: record
    description: "Batch sizing options"
    fields:
      - name: batch_size
        type: integer
        description: "The batch_size description"
        example: 1000
      - name: other_size
        type: integer
        description: "The other_size description"
        example: 2000

  - name: bwa
    type: ToolArgs
    description: "BWA aligner command line options"
    fields:
      - name: K
        type: ToolOpt
        description: "The bwa -K CLI option"
        example: 100000000
      - name: Y
        type: ToolOpt
        description: "The bwa -Y CLI option"

  - name: samtools
    type: ToolArgs
    fields:
      - name: threads
        type: ToolOpt
        ToolOpt:prefix: '-@'
        description: "Number of sort threads"
```

### Script definition

Each group is a single, flat parameter whose type is a record type — using a named `record` type or the inline destructured `record(...)` form:

```groovy
record Foo {
    batch_size: Integer = 1000
    other_size: Integer = 2000
}

params {
    output_format: String = 'bam'
    foo: Foo                                        // named record type
    bwa: ToolArgs      = record(K: ToolOpt, Y: ToolOpt)   // inline record type
    samtools: ToolArgs = record(threads: ToolOpt)
}
```

No new grammar is introduced: `params {}` is still a flat list of `name: Type = default`, and the nesting lives entirely in the type system.

### The `ToolOpt` and `ToolArgs` types

- **`ToolOpt`** — a single option; serializes per the common rules above; exposes `.name` / `.value`.
- **`ToolArgs`** — a record type whose fields are `ToolOpt`. Interpolating it renders the space-joined serialization of its **set** options, in declaration order (proposed default; see [Open Question 2](#open-questions)). Because it is an ordinary record type, fields remain accessible as `params.bwa.K`, `params.bwa.K.value`, etc.

```groovy
// params.bwa  ->  "-K 100000000 -Y"   (explicit property of the ToolArgs type)
bwa mem ${params.bwa} -t $task.cpus $index $reads
```

---

## Migration from `ext.args`

```groovy
// Legacy
process BWA_MEM {
    script:
    def args = task.ext.args ?: ''
    """ bwa mem $args -t $task.cpus $index $reads """
}
// nextflow.config
withName: 'BWA_MEM' { ext.args = '-K 100000000 -Y' }
```

```groovy
// Option A
params { bwa { K: ToolArg; Y: ToolArg } }

// Option B
params { bwa: ToolArgs = record(K: ToolOpt, Y: ToolOpt) }
```

```groovy
// nextflow.config (identical for A and B)
withName: 'BWA_MEM' {
    params.bwa.K = 100000000
    params.bwa.Y = true
}
```

In the process script, both options use `bwa mem ${params.bwa} ...`. `ext.args`/`ext.args2` are superseded by this mechanism (legacy/deprecated); they keep working during the transition, and a module may mix both while migrating.

---

## Summary & comparison

Both options deliver the same end-user surface: documented, typed, introspectable params; the same `params.<group>.<name>` override path; the same tool-option serialization. They differ in how a module author **declares** grouping.

| Dimension | Option A — `group` + `ToolArg` | Option B — record types |
|-----------|-------------------------------|-------------------------|
| New language grammar | New nested-block syntax in `params {}` | None (flat `params {}`, reuses records) |
| Concepts to learn | `group` attribute + `ToolArg` | record types + `ToolOpt`/`ToolArgs` |
| Built-in tool types | 1 (`ToolArg`) | 2 (`ToolOpt`, `ToolArgs`) |
| Simple-case verbosity | Lowest (`bwa { K: ToolArg }`) | Higher (`bwa: ToolArgs = record(K: ToolOpt)`) |
| Consistency with type system | Diverges — nesting parallels but isn't records | Aligns — one way to model structure (inputs, outputs, params) |
| Group-whole rendering (`${params.bwa}`) | Implicit rule (group-of-`ToolArg`) | Explicit type property (`ToolArgs`) |
| Relation to workflow-params ADR | Contradicts its "no nested params" non-goal | Honors it (params stay flat) |
| Reuse / composition of groups | Group is meta only — not a shareable type | Record types are named, `include`-able, duck-typed |
| `meta.yaml` shape | Flat list + `group` tag | Nested `fields:` mirroring records |
| Spec/script sync surface | Two parallel structures (YAML groups + nested block) | Type carries structure; less to drift |

**Option A pros:** lowest barrier to entry; reads naturally for tool authors; closest to the original tools-arguments draft and the nf-core mental model.
**Option A cons:** introduces nesting syntax the workflow-params ADR explicitly ruled out; creates a second way to express structure alongside records; whole-group rendering is implicit magic.

**Option B pros:** no new grammar; a single structural concept (records) across inputs/outputs/params; explicit, introspectable `ToolArgs` serialization; groups are first-class types that can be named, shared, and validated; removes `group` and its ambiguities.
**Option B cons:** higher conceptual load (records + two tool types); more verbose for the trivial single-option case; `meta.yaml` gains nested `fields:`.

**Recommendation (author's leaning):** **Option B**. It avoids new grammar, keeps a single way to model structure, and makes the tool-args concatenation an explicit type rather than an implicit rule — at the cost of asking authors to understand record types, which they will encounter anyway for typed inputs/outputs. Option A is preferable only if minimizing the learning curve for module authors is weighted above language consistency.

## Consequences

- Module parameters and tool arguments collapse into one model (`params {}`), one override path (`params.<group>.<name>`), and one spec (`meta.yaml`) — replacing `ext.args` strings and the separate `tools.*.args` namespace.
- Whichever option is chosen, the user-facing surface (CLI / config / override path / serialization) is identical, so the choice is an authoring-ergonomics and language-consistency decision, not a breaking one.
- `task.ext.args` / `ext.args2` are **superseded** by this mechanism and become legacy/deprecated; they remain functional for backward compatibility during a transition period, with a removal timeline to be decided.
- Option A introduces grammar the workflow-params ADR ruled out; Option B raises the floor of concepts an author must know (record types) but adds no new grammar.

## Open Questions

1. **Subcommand collisions** (both options): a tool with the same option name in two subcommands (`samtools view -o` vs `samtools sort -o`) cannot be distinguished within one group. Options: a group per subcommand (`params.samtools_view`), nested groups, aliasing, or scope limitation.
2. **Argument ordering** within a tool group's rendering — declaration order (proposed) vs explicit ordering.
3. **Repeated flags & conditional args** — how to model repeated flags (`-I a -I b`, possibly via a `List`-valued option) and mutually-exclusive / dependent arguments (carried over from the tools-arguments draft).
4. **Typed tool-option values** — whether to validate the underlying value beyond `enum` (e.g. `K` must be an integer).
5. **Custom type API** — finalize the `parse`/`serialize` protocol, type registration/discovery, and how `<Type>:*` spec keys are passed in and validated.
6. **`meta.yaml` ↔ script authority** — must the spec and the script `params {}` block agree exactly? Is one generated from the other, and is divergence an error?
7. **(Option B)** Guidance on inline `record(...)` vs named `record` types for parameter groups.

## Links

- Supersedes: *Module Parameters* draft — [PR #6769](https://github.com/nextflow-io/nextflow/pull/6769)
- Supersedes: *Tools Arguments* — unpublished draft (not committed to `adr/`)
- Builds on: [Workflow params](20250825-workflow-params.md), [Module System](20251114-module-system.md), [Record types](20260306-record-types.md)
