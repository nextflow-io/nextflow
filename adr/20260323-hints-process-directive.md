# `hints` process directive for executor-specific scheduling hints

- Authors: Rob Syme
- Status: draft
- Deciders: Paolo Di Tommaso, Ben Sherman, Rob Syme
- Date: 2026-03-23
- Tags: directive, executor, scheduling

Technical Story: [nextflow-io/nextflow#5917](https://github.com/nextflow-io/nextflow/issues/5917), [PR #6957](https://github.com/nextflow-io/nextflow/pull/6957)

## Summary

Introduce a `hints` process directive for executor-specific scheduling hints that don't map to existing resource directives. The first use case is AWS Batch consumable resources for license-seat-aware scheduling.

## Problem Statement

Users running commercially licensed software (e.g. DRAGEN, Schrodinger) need to limit concurrent job execution based on available license seats — not just within a single pipeline run (`maxForks`), but across multiple concurrent runs. AWS Batch introduced [resource-aware scheduling](https://docs.aws.amazon.com/batch/latest/userguide/resource-aware-scheduling.html) in February 2025, which models this natively via consumable resources. Nextflow has no way to pass this information to the executor.

More broadly, cloud batch systems expose scheduling knobs that don't fit neatly into Nextflow's existing resource directives (`cpus`, `memory`, `disk`, `accelerator`, `queue`). Today users must resort to `clusterOptions` (string-based, fragile) or custom `ext` attributes (no standard semantics). A structured, namespaced directive would provide a clean extension point.

## Goals or Decision Drivers

- Avoid a proliferation of narrow, executor-specific directives (e.g. `consumableResources`, `schedulingPolicy`, etc.)
- Provide a single extension point that executors can consume selectively
- Use a structured format (not freeform strings like `clusterOptions`)
- Keep semantics clear: hints influence scheduling/placement but do not define core resource requirements
- Distinguish from `resourceLabels` which are metadata tags, not scheduling behavior

## Non-goals

- Creating/managing cloud-side resources (e.g. AWS consumable resources) from within Nextflow
- Client-side validation of hint values against cloud APIs
- Replacing existing directives (`cpus`, `memory`, `accelerator`, `queue`) — those remain the right place for standard resources

## Considered Options

### Option 1: Dedicated `consumableResources` directive

```nextflow
process runDragen {
    consumableResources 'my-dragen-license': 1
}
```

- Good, because it's explicit and self-documenting
- Good, because the syntax maps cleanly to the AWS API
- Bad, because it adds a directive for a single executor feature
- Bad, because future executor-specific features would each need their own directive

### Option 2: Overload `resourceLabels`

```nextflow
process runDragen {
    resourceLabels 'consumable-resource:my-dragen-license': 1
}
```

- Good, because it reuses an existing directive
- Bad, because `resourceLabels` are metadata tags (write-only) — they should not influence scheduling behavior
- Bad, because it conflates two different concerns (labeling vs. scheduling)
- Bad, because `resourceLabels` uses replacement semantics in config overrides — a systems administrator setting `resourceLabels` for cost tracking in a shared compute environment would have their labels wiped out if a user sets consumable resources (and vice versa), since the entire map is replaced rather than merged

### Option 3: New `hints` directive with namespaced keys

```nextflow
process runDragen {
    hints 'consumable-resource:my-dragen-license': 1
}
```

- Good, because it provides a single, extensible directive for all executor-specific hints
- Good, because namespaced keys make intent clear and avoid collisions
- Good, because it clearly signals "this influences scheduling" without overloading existing directives
- Good, because executors can selectively consume the hints they understand and ignore the rest
- Bad, because it introduces a new directive (though one that replaces many potential future ones)

## Solution or decision outcome

Introduce a `hints` process directive (Option 3) with namespaced keys. Executors consume the hints they understand and silently ignore the rest.

## Rationale & discussion

### Syntax

The directive accepts a map of namespaced key-value pairs. It is repeatable (multiple calls accumulate) and supports config overrides:

```nextflow
// In process DSL
process runDragen {
    hints 'consumable-resource:my-dragen-license': 1
    hints 'consumable-resource:other-license': 2
    cpus 4
    memory '16 GB'

    script:
    """
    dragen --ref-dir /ref ...
    """
}

// In nextflow.config
process {
    withName: 'runDragen' {
        hints = ['consumable-resource:my-dragen-license': 1]
    }
}
```

### Namespacing

Keys use a colon-separated namespace to group related hints:

| Namespace | Key example | Value | Executor |
|-----------|-------------|-------|----------|
| `consumable-resource` | `consumable-resource:my-license` | Integer (quantity) | AWS Batch |

Future namespaces could include scheduling priorities, placement constraints, or other executor-specific features without requiring new directives.

### Executor mapping

For the initial implementation, the AWS Batch executor maps `consumable-resource:*` hints to `ConsumableResourceProperties` on `RegisterJobDefinitionRequest`:

- The portion after `consumable-resource:` becomes the resource name/ARN
- The value becomes the quantity
- Resources are set on the job definition (not as submit-time overrides)
- Different hint configurations produce distinct job definition hashes

### Important caveat: FIFO queue ordering

AWS Batch job queues use FIFO ordering by default. A job waiting for a consumable resource blocks all subsequent jobs in the same queue — even those that don't require the resource. Users should use a dedicated job queue (via the `queue` directive) or a fair-share scheduling policy.

## Links

- [AWS Batch Resource-Aware Scheduling](https://docs.aws.amazon.com/batch/latest/userguide/resource-aware-scheduling.html)
- [GitHub Issue #5917](https://github.com/nextflow-io/nextflow/issues/5917)
- [PR #6957](https://github.com/nextflow-io/nextflow/pull/6957) — initial implementation (to be reworked)
