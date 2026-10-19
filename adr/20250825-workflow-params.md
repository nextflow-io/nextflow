# Workflow params

- Authors: Ben Sherman
- Status: accepted
- Date: 2025-08-25
- Tags: lang, static-types, params

## Summary

Introduce a unified, statically typed way to declare the top-level inputs (i.e. parameters) of a workflow.

## Problem Statement

Pipeline parameters in Nextflow are currently declared using property assignments:

```groovy
params.reads = "$baseDir/data/ggal/ggal_gut_{1,2}.fq"
params.transcriptome = "$baseDir/data/ggal/ggal_1_48850000_49020000.Ggal71.500bpflank.fa"
params.multiqc = "$baseDir/multiqc"
```

This approach has several limitations:

- **No type annotations**: Parameter types cannot be expressed in the script. The type of a parameter can only be inferred from its default value, which may be ambiguous (e.g., a default value of `null`, a `String` that should be treated as a `Path`).

- **Heuristic type coercion**: When a parameter is supplied on the command line, Nextflow coerces the string value to a type using heuristics (e.g., `'true'` → boolean `true`, `'42'` → integer `42`). These heuristics are not always correct and can produce surprising values.

- **No built-in validation**: There is no built-in way to validate that a parameter is required, or that a parameter value has the correct type. The script must validate params by hand, or defer to an external JSON Schema file (`nextflow_schema.json`).

- **Scattered declarations and usage**: Parameters may be declared anywhere in the script or across multiple scripts, so there is no single view of the pipeline parameters. Parameters can also be used anywhere in the pipeline, even outside the script where they are declared, which makes it impossible to validate params usage at compile-time.

## Goals

- Declare all parameters in one place in the script, with documentation.

- Provide explicit type annotations for parameters, so that params can be validated at compile-time and supported by the IDE.

- Clearly distinguish between required and optional parameters.

- Coerce CLI parameter values based on declared types, rather than relying on heuristics.

## Non-goals

- Removing the legacy `params.foo = bar` syntax. Legacy parameters must continue to work without modification.

- Changing the `params` config scope. Params can still be declared in the config file, although some best practices apply.

- Replacing `nextflow_schema.json`. The `params` block addresses many of the same needs, but existing pipelines that use a JSON Schema should not be required to migrate. A native integration with `nextflow_schema.json` can be explored in the future.

- Supporting nested params. The `params` block only supports a flat list of params. Nested params can still be used in the config, but they do not have first-class support at this time.

## Decision

Introduce the `params` block for declaring pipeline parameters. Each parameter is declared with a name, a type, and an optional default value:

```groovy
params {
    // Path to the input samplesheet
    input: Path

    // Whether to save intermediate files
    save_intermeds: Boolean = false
}
```

Nextflow uses the declared types to validate parameter usage in the script and to coerce CLI parameter values at runtime.

## Core Capabilities

### Parameter declaration

The `params` block consists of parameter *declarations*. Each parameter is declared as `name: Type` (required) or `name: Type = default` (optional with default):

```groovy
params {
    input: Path                 // required
    extra_file: Path?           // optional (defaults to null)
    db_file: Path = 'db.json'   // optional with default
    flag: Boolean               // boolean params default to false
}
```

All standard Nextflow types except `Channel` and `Value` can be used for parameter type annotations.

### Required and optional parameters

A parameter without a default value is *required*. If a required parameter is not supplied at runtime (via the command line, a params file, or the config), the run fails immediately with an informative error.

A parameter with the `?` suffix on its type is *optional* and will be `null` if not supplied. Boolean parameters without a default value implicitly default to `false`.

### Type-based CLI coercion

When a parameter is supplied on the command line, Nextflow converts the string value to the declared type:

| Declared type | String input | Resolved value |
|---|---|---|
| `Boolean`    | `'true'`   | `true`          |
| `Integer`    | `'42'`     | `42`            |
| `Float`      | `'3.14'`   | `3.14`          |
| `Duration`   | `'1h'`     | `Duration.of('1h')` |
| `MemoryUnit` | `'8 GB'`   | `MemoryUnit.of('8 GB')` |
| `Path`       | `'/data'`  | `Path.of('/data')` |

This replaces the heuristic type detection used for legacy parameters.

### Compile-time validation

Every script in the pipeline can access legacy parameters. That flexibility comes at the cost of compile-time validation and modularity.

When a module references a param, it assumes that the workflow using it always defines that param. Nothing can check this assumption at compile-time, so a missing param fails only at runtime.

The `params` block solves this problem by defining all params in one place. It declares the inputs of the entry workflow, similar to the `take:` section in named workflows. Parameters should be passed to processes and workflows as explicit inputs, so that the type checker can validate every variable reference against a local declaration.

For example, the following workflow:

```groovy
// main.nf
params.input = '...'

workflow {
    HELLO()
}

// hello.nf
workflow HELLO {
    println "input = ${params.input}"
}
```

Can be rewritten as follows:

```groovy
// main.nf
params {
    input: String
}

workflow {
    HELLO(params.input)
}

// hello.nf
workflow HELLO {
    take:
    input: String

    main:
    println "input = ${input}"
}
```

Typed parameters can still be used globally by all scripts for backwards compatibility. However, the type checker only validates params used in the entry workflow and `output` block. To get the full benefit of type checking, pipelines should eventually be migrated as shown above.

### Script and config params

Parameters can also be defined in config files:

```groovy
params {
    outdir = 'results'
    publish_dir_mode = 'copy'
}
```

Config params continue to work as before. As a best practice, they should be used only to "configure the configuration."

Some config params can be replaced with native functionality, e.g., `outputDir` and `workflow.output.mode` for the above. The nf-core [institutional configs](https://github.com/nf-core/configs), which let users run a pipeline with their institutional config entirely from the command line, have no such replacement and are the clearest use case for config params.

Config params are also propagated to the script, because the config file can overwrite script params (e.g. in a profile). But since the script `params` block only allows params that were explicitly declared, it must distinguish config params from invalid params (e.g. a command line param with a typo).

To prevent a circular dependency between the script execution and config resolution, parameters are resolved as follows:

1. Load *CLI params* from command line, params file

2. Load config files
   - Params declared in the `params` scope are *config params*
   - If a config setting references an undeclared param, report an error
   - Params assigned in a profile are also marked as config params. They should be used to overwrite existing params or potential script params
   - CLI params override config params

3. Execute script, resolve `params` block
   - CLI params and config params override default values in `params` block
   - If a required script param is undefined, report an error
   - If a CLI param is not declared in the `params` block and is not a config param, report an error

In other words, params are applied in the following order (lowest to highest precedence):

1. Default value in the `params` block
2. Config file (`params { param = value }`)
3. Params file (`-params-file params.json`)
4. Command-line arguments (`--param value`)

Any parameter supplied via command line or params file must be declared in the script or config. Supplying an undeclared parameter is an error.

## Links

- Community issue: [#4669](https://github.com/nextflow-io/nextflow/issues/4669)
- [Workflow outputs ADR](./20251020-workflow-outputs.md)
- [Record types ADR](./20260306-record-types.md)

## Appendix

### Runtime type analysis via reflection

Validating and converting params against declared types requires the full type annotations at runtime. Parameterized types such as `List<String>` must provide both the type (`List`) and the generic type arguments (`[String]`).

During compilation, type annotations are modeled using `ClassNode`, which provides the "raw" type and type arguments via `getTypeClass() -> Class` and `getGenericsTypes() -> GenericsType[]`.

At runtime, type annotations are modeled using `Type`, which has two cases:

- If the type is parameterized, it is a `ParameterizedType`, which provides the "raw" type and type arguments via `getRawType() -> Class` and `getActualTypeArguments() -> Type[]`.

- Otherwise, the type is a `Class` corresponding to the raw type.

This type information can be obtained at runtime from the following entities:

- Class fields via `Field::getGenericType() -> Type`
- Method parameters via `Parameter::getParameterizedType() -> Type`

For this reason, the `params` block is compiled as a class, so that each parameter declaration is a field, which can model a parameterized type.

Type annotations can be marked as nullable using the `?` suffix. This marker is compiled as a custom `@Nullable` annotation on the corresponding field, so that the runtime can use this information.

For example, when loading a JSON file as a collection of records, Nextflow uses the given record type to validate each JSON object in the collection:

- String values that map to a record field with type `Path` are converted to Path values
- If a JSON object is missing a record field that is marked as nullable, it is considered valid

Every other context uses type annotations only at compile-time. Pipeline parameters need them at runtime, to validate and convert external input data to the expected type.
