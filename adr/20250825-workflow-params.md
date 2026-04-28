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

- **Heuristic type coercion**: When a parameter is supplied on the command line, Nextflow attempts to coerce the string value to the appropriate type using heuristics (e.g., `'true'` → boolean `true`, `'42'` → integer `42`). These heuristics are not always correct and can lead to unexpected behavior.

- **No built-in validation**: There is no built-in way to validate that a parameter is required, or that a parameter value has the correct type. Validation must be done manually in the script, or through an external JSON Schema file (`nextflow_schema.json`).

- **Scattered declarations and usage**: Parameters may be declared anywhere in the script or across multiple scripts, making it difficult to get a single view of the pipeline parameters. Parameters can be used anywhere in the pipeline, even outside the script where they are declared, making it impossible to validate params usage at compile-time.

## Goals

- Declare all parameters in one place in the script, with documentation.

- Provide explicit type annotations for parameters, enabling compile-time validation and IDE support.

- Clearly distinguish between required and optional parameters.

- Coerce CLI parameter values based on declared types, rather than relying on heuristics.

- Support collection-type parameters that can be loaded from structured files (CSV, JSON, YAML).

## Non-goals

- Removing the legacy `params.foo = bar` syntax -- legacy parameters must continue to work without modification.

- Changing the `params` config scope -- params can still be declared in the config file, although some best practices apply.

- Replacing `nextflow_schema.json` -- while the `params` block addresses many of the same needs, existing pipelines that use a JSON Schema should not be required to migrate. A native integration with `nextflow_schema.json` can be explored in the future.

- Supporting nested params -- the `params` block only supports a flat list of params. Nested params can still be used in the config, but they do not have first-class support at this time.

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

Typed parameters are used to validate parameter usage in the script, and to coerce CLI parameter values at runtime.

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

### Samplesheets as collection-type parameters

A parameter with a collection type (`List`, `Set`, `Bag`) can be supplied as a file path. Nextflow parses the file and assigns the resulting collection to the parameter. Supported formats are CSV, JSON, and YAML:

```groovy
params {
    samples: List<Sample>   // can be supplied as a CSV, JSON, or YAML file path
}

record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path
}
```

The file contents must be compatible with the declared element type; an error is thrown if they are not. CSV files must include a header row and use a comma as the column separator.

The collection-type parameter can use a generic type such as `Map` or `Record`, or a custom record type to enable further validation. In the above example, using the `Sample` type ensures that each samplesheet row is validated against the record fields and the `fastq_1` and `fastq_2` columns are treated as file paths.

This feature allows collection-type parameters to serve as *samplesheet inputs*, which simplifies the workflow logic and allows it to be agnostic to the input format:

```groovy
// before (CSV only)
ch_samples = channel.fromPath(param.samples)
    .flatMap { csv ->
        csv.splitCsv(header: true, sep: ',')
    }
    .map { r ->
        record(id: r.id, fastq_1: file(r.fastq_1), fastq_2: file(r.fastq_2))
    }

// after (CSV, JSON, or YAML)
ch_samples = channel.fromList(param.samples)
```

### Compile-time validation

Legacy parameters can be accessed globally by all scripts in the pipeline. While this approach is flexible, it prevents compile-time validation and breaks modularity.

When a module references a param, it implicitly assumes that the param will always be defined by the workflow that uses it. This assumption cannot be validated at compile-time, so if the param is missing, an error will occur only at runtime.

The `params` block solves this problem by defining all params in one place. It serves as the inputs for the entry workflow, similar to the `take:` section in named workflows. Parameters should be passed to processes and workflows as explicit inputs, so that every variable reference can be validated against local declarations.

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

Typed parameters can still be used globally by all scripts for backwards compatibility. However, the type checker will only validate params used in the entry workflow and `output` block. Users should eventually migrate their pipelines as shown above for effective type checking.

### Script and config params

Parameters can also be defined in config files:

```groovy
params {
    outdir = 'results'
    publish_dir_mode = 'copy'
}
```

Config params continue to work as before. As a best practice, they should be used only to "configure the configuration."

Some config params can be replaced with native functionality, e.g., `outputDir` and `workflow.output.mode` for the above. The nf-core [institutional configs](https://github.com/nf-core/configs), which enable users to run a pipeline with their institutional config entirely from the command line, cannot be easily replaced and provide a clear use case for config params.

Config params are also propagated to the script since the config file can overwrite script params (e.g. in a profile). However, since the script `params` block only allows params that were explicitly declared, it needs to be able to distinguish between config params and invalid params (e.g. command line param with a typo).

To prevent a circular dependency between the script execution and config resolution, parameters are resolved as follows:

1. Load *CLI params* from command line, params file

2. Load config files
   - Params declared in the `params` scope are *config params*
   - If a config setting references an undeclared param, report an error
   - Params assigned in a profile are also marked as config params -- they should be used to overwrite existing params or potential script params
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

Validating and converting params against declared types requires the type annotations to be fully available at runtime. Parameterized types such as `List<String>` must provide both the type (`List`) and the generic type arguments (`[String]`).

During compilation, type annotations are modeled using `ClassNode`, which provides the "raw" type and type arguments via `getTypeClass() -> Class` and `getGenericsTypes() -> GenericsType[]`.

At runtime, type annotations are modeled using `Type`, for which there are two primary cases:

- If the type is parameterized, it is a `ParameterizedType`, which provides the "raw" type and type arguments via `getRawType() -> Class` and `getActualTypeArguments() -> Type[]`.

- Otherwise, the type is a `Class` corresponding to the raw type.

This type information can be obtained at runtime from the following entities:

- Class fields via `Field::getGenericType() -> Type`
- Method parameters via `Parameter::getParameterizedType() -> Type`

For this reason, the `params` block is compiled as a class, so that each parameter declaration is a field which can model a parameterized type.

Type annotations can be marked as nullable using the `?` suffix. This marker is compiled as a custom `@Nullable` annotation on the corresponding field, so that the runtime can use this information.

For example, when loading a JSON file as a collection of records, Nextflow uses the given record type to validate each JSON object in the collection:

- String values that map to a record field with type `Path` are converted to Path values
- If a JSON object is missing a record field that is marked as nullable, it is considered valid

While type annotations are used only at compile-time in all other contexts, they are needed at runtime for pipeline parameters in order to validate and convert external input data to the expected type.

### Standard library functions for loading structured files

The automatic loading of samplesheet inputs is supported only in the `params` block. It could also be useful to load structured files with functions. For example:

```groovy
samples = fromJson('samples.json')
```

Where `fromJson` is a function that loads arbitrary data from a JSON file.

A function is more flexible because it can be used anywhere in pipeline code, whereas the automatic samplesheet loading can only be used for pipeline-level inputs. For example, a process might produce a JSON file that needs to be read in workflow logic and processed by downstream tasks in parallel.

However, data-loading functions like `fromJson` cannot be statically typed, since the data file could contain anything (number, string, list, map, etc). Primitive values can be coerced using a Groovy-style cast (e.g. `fromJson('...') as Map`), but this approach does not support parameterized types or Nextflow record types. A function like `fromJson` also assumes a specific file format, which is overly restrictive for a pipeline-level input that could be supplied in a variety of formats.

Loading structured files with a typed parameter such as `samples: List<Sample>` allows the input to have a well-defined type at compile-time and allows it to be sourced from any data format. While currently only CSV, JSON, and YAML are supported, the format can be made extensible in the future so that users can integrate their own data formats (e.g. Parquet) via plugins. Data-loading functions like `fromJson` can also be implemented for other use cases in the future, but they are not sufficient for handling pipeline-level inputs.

### Typed parameters and schemas

While this ADR does not specify any native integration with JSON schema, it is worth addressing how typed parameters are expected to interact with schemas in the future.

A common approach for Nextflow pipelines is to define a JSON schema for the pipeline parameters. A samplesheet param is typically defined as a file input with its own *samplesheet schema*. Developers typically load and validate samplesheet params using the `samplesheetToList` function from the `nf-schema` plugin:

```groovy
include { samplesheetToList } from 'plugin/nf-schema'

params.input = null

workflow {
    samples = samplesheetToList(params.input, "assets/schema_input.json")
    channel.fromList(samples).view()
}
```

The parameter schema provides similar validation and type coercion as a typed parameter. For example, the following record type:

```groovy
record Sample {
    id: String
    fastq_1: Path
    fastq_2: Path?
}
```

Is equivalent to the following JSON schema:

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "type": "array",
  "items": {
    "type": "object",
    "properties": {
      "id": { "type": "string" },
      "fastq_1": { "type": "string", "format": "file-path", "exists": true },
      "fastq_2": { "type": "string", "format": "file-path", "exists": true }
    },
    "required": ["id", "fastq_1"]
  }
}
```

The `samplesheetToList` function has the same limitation as `fromJson` described above -- it is not statically typed. Even with the schema file, the type checker cannot guarantee a specific return type at compile-time because there is no type annotation. The above example must be rewritten as follows in order to be statically typed:

```groovy
params {
    input: List<Sample>
}

workflow {
    channel.fromList(params.input).view()
}
```

On the other hand, the samplesheet schema can specify additional validations that cannot be expressed with Nextflow types, such as min/max constraints for numbers and pattern constraints for strings.

It is not yet clear whether it is better for Nextflow to enforce these schema properties at runtime, or for users to implement the equivalent validation in their pipeline code. If needed, Nextflow should be able to augment the `params` block with the parameter schema to provide this extra validation. The above example would work without any additional code changes.
