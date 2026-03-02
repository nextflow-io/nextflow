# Tools Arguments

- Authors: Paolo Di Tommaso
- Status: draft
- Date: 2026-01-28
- Tags: modules, tools, arguments, configuration

## Context and Problem Statement

The current adoption of `task.ext.args` for passing tool arguments in Nextflow modules is extremely convoluted. This pattern makes it impossible to have a consistent, programmatic definition for parameters exposed by a module.

**Current approach (problematic):**
```groovy
process BWA_MEM {
    script:
    def args = task.ext.args ?: ''
    """
    bwa mem $args -t $task.cpus $index $reads
    """
}
```

```groovy
// Configuration
withName: 'BWA_MEM' {
    ext.args = '-K 100000000 -Y'
}
```

**Problems:**
- Arguments are opaque strings with no validation
- No documentation of available options
- No type safety or IDE support
- Easy to introduce typos or invalid combinations
- Impossible to programmatically introspect module capabilities

## Decision

Extend the `tools` definition in the module spec (`meta.yaml`) to support an `args` component that declares exposed tool command line arguments.

**Key requirement:** The argument name must match the tool's actual option name.

## Tool Arguments Specification

### Definition in meta.yaml

```yaml
tools:
  - samtools:
      description: SAMtools
      homepage: http://www.htslib.org/
      args:
        output_fmt:
          type: string
          enum: ["sam", "bam", "cram"]
          description: "Output format"

  - bwa:
      description: BWA aligner
      homepage: http://bio-bwa.sourceforge.net/
      args:
        K:
          type: integer
          description: "Process INT input bases in each batch"
          prefix: '-'
        Y:
          type: boolean
          description: "Use soft clipping for supplementary alignments"
          prefix: '-'
```

### Argument Attributes

| Attribute | Required | Default | Description |
|-----------|----------|---------|-------------|
| `type` | No | `string` | Data type: `boolean`, `integer`, `float`, `string` |
| `description` | No | - | Human-readable description |
| `enum` | No | - | List of allowed values |
| `prefix` | No | `'--'` | Command-line prefix for the argument |
| `default` | No | - | Default value if not configured |

### Prefix Resolution

The argument prefix determines how the argument is formatted on the command line:

| Argument | Prefix | Value | Resolved |
|----------|--------|-------|----------|
| `output_fmt` | `--` (default) | `"cram"` | `--output_fmt cram` |
| `K` | `-` | `100000000` | `-K 100000000` |
| `Y` | `-` | `true` | `-Y` |

**Boolean arguments:** When `type: boolean` and value is `true`, only the prefix+name is emitted (no value).

## Script Usage

In module scripts, access arguments via the `tools` implicit variable:

```groovy
process BWA_MEM {
    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.bam"), emit: aligned

    script:
    // tools.bwa.args.K → "-K 100000000"
    // tools.bwa.args.Y → "-Y"
    // tools.bwa.args   → "-K 100000000 -Y"  (all args concatenated)
    """
    bwa mem ${tools.bwa.args} -t $task.cpus $index $reads \
        | samtools sort ${tools.samtools.args} -o ${prefix}.bam -
    """
}
```

### Implicit Variable Reference

| Expression | Description | Example Output |
|------------|-------------|----------------|
| `tools.<tool>.args.<arg>` | Single formatted argument | `"-K 100000000"` |
| `tools.<tool>.args` | All arguments concatenated | `"-K 100000000 -Y"` |

## Configuration Usage

Arguments are configured using `tools.<toolname>.args.<argname>`:

```groovy
process {
    withName: 'BWA_MEM' {
        tools.bwa.args.K = 100000000
        tools.bwa.args.Y = true
        tools.samtools.args.output_fmt = "cram"
    }
}
```

## CLI Usage

Arguments can be overridden via CLI:

```bash
# For module run command
nextflow module run nf-core/bwa-align \
    --tools.bwa.K=100000000 \
    --tools.bwa.Y \
    --tools.samtools.output_fmt=cram

# For standard workflow execution
nextflow run <script> \
    -process.BWA_MEM.tools.bwa.K=100000000 \
    -process.BWA_MEM.tools.bwa.Y \
    -process.BWA_MEM.tools.samtools.output_fmt=cram
```

## Benefits

| Aspect | `ext.args` (Legacy) | Tool Arguments (New) |
|--------|---------------------|----------------------|
| Documentation | None | In meta.yaml |
| Type Safety | None | Validated |
| IDE Support | None | Autocompletion |
| Introspection | Impossible | Programmatic access |
| Validation | None | Type + enum validation |
| Clarity | Opaque strings | Named, typed arguments |

## Complete Example

### meta.yaml

```yaml
name: nf-core/bwa-align
version: "1.2.4"
description: Align reads to reference genome using BWA-MEM algorithm

requires:
  nextflow: ">=24.04.0"

tools:
  - bwa:
      description: |
        BWA is a software package for mapping DNA sequences
        against a large reference genome.
      homepage: http://bio-bwa.sourceforge.net/
      documentation: https://bio-bwa.sourceforge.net/bwa.shtml
      doi: 10.1093/bioinformatics/btp324
      license: ["GPL-3.0-or-later"]
      identifier: biotools:bwa
      args:
        K:
          type: integer
          description: "Process INT input bases in each batch"
          prefix: '-'
        Y:
          type: boolean
          description: "Use soft clipping for supplementary alignments"
          prefix: '-'
        T:
          type: integer
          description: "Minimum score to output"
          prefix: '-'
          default: 30

  - samtools:
      description: Tools for manipulating alignments in SAM/BAM format
      homepage: http://www.htslib.org/
      license: ["MIT"]
      identifier: biotools:samtools
      args:
        output_fmt:
          type: string
          enum: ["sam", "bam", "cram"]
          description: "Output format"
          default: "bam"
        threads:
          type: integer
          description: "Number of threads"
          prefix: '-@'
```

### main.nf

```groovy
process BWA_MEM {
    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(index)

    output:
    tuple val(meta), path("*.bam"), emit: aligned
    path "versions.yml", emit: versions

    script:
    def prefix = meta.id
    """
    bwa mem ${tools.bwa.args} -t $task.cpus $index $reads \
        | samtools sort ${tools.samtools.args} -o ${prefix}.bam -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(bwa 2>&1 | grep Version | cut -d' ' -f2)
        samtools: \$(samtools --version | head -1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
```

### Configuration

```groovy
process {
    withName: 'BWA_MEM' {
        tools.bwa.args.K = 100000000
        tools.bwa.args.Y = true
        tools.samtools.args.output_fmt = "cram"
    }
}
```

## Migration from ext.args

### Before (Legacy)

```groovy
// main.nf
process BWA_MEM {
    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    """
    bwa mem $args -t $task.cpus $index $reads | samtools sort $args2 -o out.bam -
    """
}

// nextflow.config
withName: 'BWA_MEM' {
    ext.args = '-K 100000000 -Y'
    ext.args2 = '--output-fmt cram'
}
```

### After (Tool Arguments)

```groovy
// main.nf
process BWA_MEM {
    script:
    """
    bwa mem ${tools.bwa.args} -t $task.cpus $index $reads \
        | samtools sort ${tools.samtools.args} -o out.bam -
    """
}

// nextflow.config
withName: 'BWA_MEM' {
    tools.bwa.args.K = 100000000
    tools.bwa.args.Y = true
    tools.samtools.args.output_fmt = "cram"
}
```

## Rationale

**Why argument names must match tool options?**
- Eliminates mapping layer between config and actual tool
- Users familiar with the tool can immediately understand the config
- Documentation can link directly to tool manuals
- Reduces cognitive overhead and potential for errors

**Why per-tool namespacing?**
- Modules often wrap multiple tools (e.g., bwa + samtools)
- Clear separation prevents argument name collisions
- Enables tool-specific documentation and validation

**Why typed arguments?**
- Catch errors at configuration time, not runtime
- Enable IDE autocompletion and validation
- Support for enums provides constrained choices

## Open Questions

1. **Argument ordering**: Should the order of `tools.<tool>.args` output match declaration order in meta.yaml or be alphabetical?

2. **Conditional arguments**: How to handle arguments that are mutually exclusive or have dependencies?

3. **Complex argument formats**: How to handle arguments with complex formats like `-I file1 -I file2` (repeated flags)?

## Open Problems

### Subcommand argument collision

A tool having the same option name in two different subcommands cannot be managed with the current design. Arguments are defined at the tool level, not at the subcommand level.

**Example:** `samtools view` and `samtools sort` both have a `-o` option with different semantics:

```bash
samtools view -o output.bam input.bam    # output file for view
samtools sort -o output.bam input.bam    # output file for sort
```

With the current spec, there's no way to distinguish between these:

```yaml
tools:
  - samtools:
      args:
        o:
          type: string
          prefix: '-'
          description: "Output file"  # Which subcommand?
```

**Potential solutions (not addressed in this ADR):**
- Subcommand namespacing: `tools.samtools.view.args.o` vs `tools.samtools.sort.args.o`
- Argument aliasing: Define tool-level unique names that map to subcommand options
- Scope limitation: Only support arguments that are consistent across subcommands
