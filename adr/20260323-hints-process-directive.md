# `hints` process directive for executor-specific scheduling hints

- Authors: Rob Syme
- Status: accepted
- Deciders: Paolo Di Tommaso, Ben Sherman, Rob Syme
- Date: 2026-03-23
- Tags: directive, executor, scheduling

## Summary

Introduce a `hints` process directive for executor-specific scheduling hints that don't map to existing directives.

## Problem Statement

Many executors can be configured in various ways on a per-task basis. For example:

- AWS Batch jobs can use *consumable resources* to limit concurrent job execution based on non-standard resources such as software license seats.

- Google Batch jobs can specify a *provisioning model* to control the use of spot vs on-demand VMs on a per-task basis.

- Seqera Scheduler supports a variety of resource and scheduling settings, including spot/on-demand provisioning.

These settings can be exposed by Nextflow as executor-specific config options, such as `google.batch.spot`, but config options are applied globally. In order to apply a setting to specific processes or tasks, it must be exposed as a process directive.

Process directives in Nextflow aim to provide a common vocabulary for executing tasks in many different environments. Directives such as `cpus`, `memory`, and `time` have broadly the same meaning across most executors, making it easier for users to write portable pipelines.

At the same time, many executors have custom settings not shared by other executors, and it is not practical to create a new process directive for every new setting. There are over 40 [process directives](https://docs.seqera.io/nextflow/reference/process#directives) at the time of writing, and every new directive adds cognitive load when a user is trying to find the right directive for a given situation.

There exist a few generic process directives already:

- The `clusterOptions` directive can be used to specify command-line arguments, primarily for HPC schedulers
- The `ext` directive supports arbitrary key-values, but is designed primarily to customize the task script (e.g. tool arguments), not executor behavior
- The `resourceLabels` directive also supports arbitrary key-values, but is intended for tagging and tracking resources, not controlling them

A new directive is needed to support executor-specific settings at a per-task level in a structured manner, without bloating the process directives for every new custom setting.

## Goals

- Provide a way to apply executor-specific settings to individual processes or tasks

- Avoid the proliferation of narrow, executor-specific directives (e.g. `consumableResources`, `schedulingPolicy`, etc.)

- Provide a single extension point that executors can consume selectively

- Allow settings to be specified as key-values, providing validation where possible

## Non-goals

- Replacing existing directives (`cpus`, `memory`, `accelerator`, `queue`) — those remain the right place for standard resources

## Decision

Introduce a `hints` process directive with namespaced keys. Executors consume the hints they understand and silently ignore the rest.

## Core Capabilities

### Syntax

The `hints` directive accepts a map of key-value pairs:

```groovy
// process definition
process runDragen {
    cpus 4
    memory '16 GB'
    hints consumableResources: ['my-dragen-license': 1, 'other-license': 2]

    script:
    """
    dragen --ref-dir /ref ...
    """
}
```

```groovy
// process config
process {
    withName: 'runDragen' {
        hints = [
            consumableResources: ['my-dragen-license': 1, 'other-license': 2]
        ]
    }
}
```

Keys are strings. Values may be any raw data type: strings, numbers, booleans, lists, or maps. Executors are responsible for defining which hints they recognize and what value type each hint expects.

In the above example, the `consumableResources` hint is given as a map of resource name to quantity. The AWS Batch executor supplies it to each job request using `ConsumableResourceProperties`.

### Namespacing

Keys can use dot-separated scopes to namespace settings as needed:

```groovy
hints consumableResources: ['my-dragen-license': 1]
hints 'scheduling.priority': 10
hints 'scheduling.provisioningModel': 'spot'
```

Keys can be routed to specific executors by prefixing with the executor name and a slash (`/`):

```groovy
hints 'awsbatch/consumableResources': ['my-dragen-license': 1]
hints 'seqera/scheduling.provisioningModel': 'spot'
hints 'k8s/nodeSelector': 'gpu=true'
```

The executor prefix gives pipeline developers the ability to target specific executors and have assurance that it won't accidentally apply to other executors (e.g. if another executor adds support for the same hint in the future).

### Validation

Nextflow should validate hints to the best of its ability, to catch errors such as typos:

- **Prefixed hints** can be validated against the set of hints declared by the corresponding executor. Unrecognized hints should be reported as errors.

- **Unprefixed hints** can be validated against the union of hints declared by all executors. Since unprefixed hints might be supported by executors that aren't currently loaded, unrecognized hints should be reported as warnings.

### Multiple hint resolution

The `hints` directive uses *replacement semantics* when specified multiple times, meaning that each `hints` setting completely replaces any previous settings:

```groovy
process {
    // generic hint
    hints = [provisioningModel: 'spot']

    // specific hint replaces generic hint
    withLabel: 'dragen' {
        hints = [consumableResources: ['my-dragen-license': 1]]
    }
}
```

Within a process definition, the `hints` directive uses *accumulation semantics*, meaning that subsequent `hints` directives are accumulated:

```groovy
process runDragen {
    // multiple separate hints
    hints provisioningModel: 'spot'
    hints consumableResources: ['my-dragen-license': 1, 'other-license': 2]

    // equivalent to...
    hints (
        provisioningModel: 'spot',
        consumableResources: ['my-dragen-license': 1, 'other-license': 2]
    )

    // ...
}
```

This behavior is consistent with other directives such as `pod` and `resourceLabels`. In practice, this means that a given `hints` setting should specify all relevant hints for the given context.

For example, the `withLabel` selector above should also specify the `provisioningModel` hint if the intention is to preserve that hint for the selected processes:

```groovy
process {
    hints = [provisioningModel: 'spot']

    withLabel: 'dragen' {
        hints = [provisioningModel: 'spot', consumableResources: ['my-dragen-license': 1]]
    }
}
```

While this approach may lead to duplication, it gives users and developers more control over which hints are applied in a given context.

### Initial hint catalog

The following hints should be supported initially:

| Hint name | Value type | Executors | Use case |
|--|--|--|--|
| `consumableResources` | `Map<String, Integer>` | AWS Batch | License-aware scheduling ([#5917](https://github.com/nextflow-io/nextflow/issues/5917)) |
| `scheduling.priority` | `Integer` | AWS Batch | Job scheduling priority ([#6998](https://github.com/nextflow-io/nextflow/issues/6998)) |
| `scheduling.provisioningModel` | `String` | Google Batch | Spot VM scheduling ([#3530](https://github.com/nextflow-io/nextflow/issues/3530)) |

## Links

- [Community issue](https://github.com/nextflow-io/nextflow/issues/5917)
