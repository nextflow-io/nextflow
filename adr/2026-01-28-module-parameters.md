# Module Parameters

- Authors: Paolo Di Tommaso
- Status: draft
- Date: 2025-01-28
- Tags: modules, parameters, configuration
- Related: [Module System ADR](20251114-module-system.md)

## Context

The module system introduces a structured approach to module configuration through parameters defined in `meta.yaml`. Parameters provide a documented, type-safe way to customize module behavior.

This document describes two aspects of module parameters:
1. **Specification** (`meta.yaml`) - declares parameter metadata for documentation, validation, and IDE support
2. **Definition in script** - how parameters are declared and used in the Nextflow module script

## Parameter Specification (meta.yaml)

Modules declare available parameters in `meta.yaml` under the `params` section. Each parameter has a name and optional attributes for type and description.

```yaml
params:
  - name: batch_size
    type: integer
    description: "Process INT input bases in each batch"
    example: 100000000

  - name: use_soft_clipping
    type: boolean
    description: "Use soft clipping for supplementary alignments"

  - name: output_format
    type: string
    description: "Output format (sam, bam, or cram)"
```

**Specification Attributes:**

| Attribute | Required | Description |
|-----------|----------|-------------|
| `name` | Yes | Parameter identifier |
| `type` | No | Data type: `boolean`, `integer`, `float`, `string`, `file`, `path` |
| `description` | No | Human-readable description |
| `example` | No | Example value for the parameter |

## Parameter Definition in Script

Parameters declared in `meta.yaml` must also be defined in the module script. Two options are considered:

### Option A: Process-level params block

Parameters are defined within the process definition using a new `params:` block:

```groovy
process BWA_MEM {
    params:
    batch_size: Integer
    use_soft_clipping: Boolean = false
    output_format: String = 'bam'

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.${params.output_format}"), emit: aligned

    script:
    def batch_arg = params.batch_size ? "-K ${params.batch_size}" : ''
    def soft_clip = params.use_soft_clipping ? "-Y" : ''
    """
    bwa mem ${batch_arg} ${soft_clip} -t $task.cpus $index $reads \
        | samtools sort --output-fmt ${params.output_format} -o ${meta.id}.${params.output_format} -
    """
}
```

**Note:** Module-scoped `params {}` blocks are still allowed and represent workflow/global parameters (same value shared across all modules).

**Pros:**
- Parameters are scoped to the process - clear ownership
- Natural extension of process definition syntax
- Unambiguous naming for config and CLI overrides
- Clear distinction between process params and global params

**Cons:**
- New syntax block in process definition
- Parameters not accessible outside the process

---

### Option B: Module-level params block

Parameters are defined at the module script level using a `params {}` block:

```groovy
params {
    reads: Path
    index: Path
    batch_size: Integer
    use_soft_clipping: Boolean = false
    output_format: String = 'bam'
}

process BWA_MEM {
    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("*.${params.output_format}"), emit: aligned

    script:
    def batch_arg = params.batch_size ? "-K ${params.batch_size}" : ''
    def soft_clip = params.use_soft_clipping ? "-Y" : ''
    """
    bwa mem ${batch_arg} ${soft_clip} -t $task.cpus $index $reads \
        | samtools sort --output-fmt ${params.output_format} -o ${meta.id}.${params.output_format} -
    """
}

workflow {
    BWA_MEM(params.reads, params.index)
}
```

**Pros:**
- Consistent with current Nextflow params syntax
- Parameters accessible in entry workflow
- Familiar pattern for existing users

**Cons:**
- Ambiguous naming when overriding from CLI and config (module name vs process name?)
- Parameters not clearly scoped to a specific process

---

## Parameter Semantics

Regardless of the definition option chosen:

**Typing and Defaults:**
- Parameters can be typed (`String`, `Integer`, `Boolean`, `Float`, `Path`, `File`)
- Default values must be literals, constants, or references to other parameters
- Type validation occurs at resolution time

**Parameter vs Input:**
| Aspect | Parameter | Input |
|--------|-----------|-------|
| Resolution | Before process execution | Per task execution |
| Mutability | Immutable once resolved | Different value per task |
| Source | Config, CLI, defaults | Channels, values |

**Usage:**
- In script context (like inputs)
- In input definitions (to define default values)

## Configuration Override

Parameters can be overridden via config:

```groovy
process {
    withName: 'BWA_MEM' {
        params.batch_size = 100000000
        params.use_soft_clipping = true
        params.output_format = 'cram'
    }
}
```

## CLI Override

Parameters can be overridden via CLI:

```bash
# For standard workflow execution
nextflow run <script> -process.BWA_MEM.params.batch_size=100000000

# For module run command
nextflow module run nf-core/bwa-align --batch_size=100000000 --use_soft_clipping
```

## Benefits

| Aspect | `ext.args` (Legacy) | Module Parameters (New) |
|--------|---------------------|-------------------------|
| Documentation | None | In meta.yaml |
| Type Safety | None | Validated |
| IDE Support | None | Autocompletion |
| Clarity | Opaque strings | Named parameters |
| Defaults | Manual | Schema-defined |

## Complete Example

```yaml
name: nf-core/bwa-align
version: "1.2.4"
description: Align reads to reference genome using BWA-MEM algorithm

requires:
  nextflow: ">=24.04.0"

params:
  - name: batch_size
    type: integer
    description: "Process INT input bases in each batch"
    example: 100000000

  - name: use_soft_clipping
    type: boolean
    description: "Use soft clipping for supplementary alignments"

  - name: output_format
    type: string
    description: "Output format (sam, bam, or cram)"

tools:
  - bwa:
      description: BWA aligner
      homepage: http://bio-bwa.sourceforge.net/
      license: ["GPL-3.0-or-later"]
      identifier: biotools:bwa
```
