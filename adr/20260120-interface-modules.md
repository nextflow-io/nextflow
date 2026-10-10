# Interface Modules for Runtime Implementation Selection

- Authors: Paolo Di Tommaso
- Status: draft
- Date: 2026-01-20
- Tags: modules, interfaces, polymorphism, benchmarking

Technical Story: https://github.com/nextflow-io/schemas/issues/11

## Summary

Introduce interface modules - abstract module definitions that declare input/output contracts without implementation. Concrete modules can implement these interfaces, and users select the target implementation at runtime via configuration. This enables tool benchmarking, user customization, and LLM-assisted pipeline construction.

## Problem Statement

The nf-core community has identified patterns where multiple tools perform equivalent operations with compatible input/output signatures. For example, multiple sequence alignment tools (ClustalO, FAMSA, KALIGN, etc.) all accept FASTA sequences and produce alignment files.

Currently, supporting runtime tool selection requires verbose boilerplate using branch operators to route data to different implementations and mix outputs back together. This approach requires updating the subworkflow for every new tool, provides no validation that tools conform to expected contracts, and is difficult for LLMs to reason about.

## Goals

- **Minimal syntax impact**: No new DSL keywords; leverage existing module system and `meta.yaml` schema
- **Familiar patterns**: Interface modules are included and used like regular modules
- **Configuration-based binding**: Implementation selection via `nextflow.config`, not workflow code
- **Contract validation**: Nextflow validates that implementations match interface contracts

## Non-goals

- Automatic implementation discovery via runtime registry queries
- Multiple inheritance or interface inheritance
- Automatic output format normalization between implementations

## Solution

### Interface Module Definition

An interface module has `type: interface` in its `meta.yaml` and no script implementation:

```yaml
# modules/@nf-core/msa-alignment/meta.yaml
name: nf-core/msa-alignment
version: "1.0.0"
description: Perform multiple sequence alignment
type: interface

keywords:
  - alignment
  - msa

ontology:
  operation: "http://edamontology.org/operation_0492"

input:
  - - meta:
        type: map
        description: Sample metadata
    - fasta:
        type: file
        description: Input sequences in FASTA format
        pattern: "*.{fa,fasta}"

output:
  alignment:
    - - meta:
          type: map
      - "*.aln.gz":
          type: file
          description: Alignment file
          pattern: "*.aln.gz"
  versions:
    - versions.yml:
        type: file
```

An optional `main.nf` stub can provide IDE support:

```groovy
process MSA_ALIGNMENT {
    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.aln.gz"), emit: alignment
    path "versions.yml",              emit: versions
}
```

### Implementation Declaration

Concrete modules declare interface compliance via `implements` in `meta.yaml`:

```yaml
# modules/@nf-core/clustalo-align/meta.yaml
name: nf-core/clustalo-align
version: "1.1.0"
description: Multiple sequence alignment using Clustal Omega

implements: nf-core/msa-alignment@>=1.0.0

tools:
  - clustalo:
      description: Clustal Omega
      homepage: http://www.clustal.org/omega/

input:
  - - meta:
        type: map
    - fasta:
        type: file
        pattern: "*.{fa,fasta}"

output:
  alignment:
    - - meta:
          type: map
      - "*.aln.gz":
          type: file
  versions:
    - versions.yml:
        type: file
```

**Validation rules:**
- All interface inputs must be present in the implementation
- All interface outputs must be present with matching emit names
- Additional inputs/outputs are allowed (interface is a minimum contract)

### Using Interface Modules

Interface modules are included using standard syntax:

```groovy
include { MSA_ALIGNMENT } from '@nf-core/msa-alignment'

workflow {
    ch_fasta = Channel.fromPath(params.input)
        .map { file -> [[id: file.baseName], file] }

    MSA_ALIGNMENT(ch_fasta)

    MSA_ALIGNMENT.out.alignment.view()
}
```

### Implementation Resolution

Implementation binding is configured in `nextflow.config`:

```groovy
modules {
    '@nf-core/msa-alignment' = '1.0.0'
    '@nf-core/clustalo-align' = '1.1.0'
    '@nf-core/famsa-align' = '2.0.0'

    interfaces {
        'nf-core/msa-alignment' = 'nf-core/clustalo-align'
    }
}
```

**Resolution strategies:**

1. **Static binding**:
   ```groovy
   interfaces {
       'nf-core/msa-alignment' = 'nf-core/clustalo-align'
   }
   ```

2. **Parameter-based binding**:
   ```groovy
   interfaces {
       'nf-core/msa-alignment' = params.aligner
   }
   ```

3. **Per-sample binding** (closure):
   ```groovy
   interfaces {
       'nf-core/msa-alignment' {
           resolve = { meta, fasta -> meta.tool }
       }
   }
   ```

### Runtime Behavior

When Nextflow encounters an interface module invocation:

1. **Parse time**: Load interface module, identify `type: interface`
2. **Resolution**: Look up bound implementation in `modules.interfaces`
3. **Validation**: Verify implementation declares `implements: <interface>`
4. **Dispatch**: Route invocation to the concrete implementation

For per-sample resolution, each channel item is evaluated against the resolve closure, routed to appropriate implementations, and outputs are collected into unified channels.

### Schema Extensions

Add to module spec schema:

| Field | Type | Description |
|-------|------|-------------|
| `type` | `enum["process", "workflow", "interface"]` | Module type; `interface` for abstract contracts |
| `implements` | `string` | Interface module this implements (with optional version constraint) |
| `ontology.operation` | `uri` | EDAM operation URI for semantic discovery |

### CLI Extensions

```bash
nextflow module list --implements nf-core/msa-alignment
nextflow module search --type interface
nextflow module validate nf-core/clustalo-align
```

## Links

- Related: [Module System ADR](20251114-module-system.md)
- Discussion: [nextflow-io/schemas#11](https://github.com/nextflow-io/schemas/issues/11)
- Reference: [nf-core class-modules](https://github.com/mirpedrol/class-modules)

## Open Questions

1. **Default implementation**: Should interfaces specify a default if none is configured?
2. **Multiple interfaces**: Can a module implement more than one interface?
3. **Interface versioning**: What constitutes a breaking change to an interface?
