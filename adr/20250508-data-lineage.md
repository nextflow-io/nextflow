# Data lineage

- Authors: Ben Sherman
- Status: accepted
- Deciders: Paolo Di Tommaso, Jorge Ejarque, Ben Sherman
- Date: 2025-05-08
- Tags: lineage, provenance

## Summary

Introduce native provenance tracking for files produced by Nextflow pipelines.

## Problem Statement

Scientific workflows are run many times across different environments, parameter sets, and software versions. Without provenance tracking, it is difficult to answer questions such as:

- Which exact version of a pipeline and which inputs produced a given result file?
- Why was a previously cached task re-executed on a resumed run?
- What are all the files associated with a particular sample ID or label?

Nextflow already produces this provenance information at runtime, saving some of it to the metadata cache and the work directory, but much of this information is discarded at the end of a run. Even the cache and work directory are ephemeral by design, not suitable as a persistent store of lineage data.

Workflow outputs are typically published to an external output directory, organized by embedding metadata into file paths (e.g. a sub-directory for each sample ID). This approach is arbitrary and error-prone, as it relies on users already knowing how outputs are organized for a given run. The `output` block with index files makes this problem more manageable, but does not completely solve it.

## Goals

- Record every workflow run, task run, and output file produced by Nextflow into a lineage graph.

- Define an immutable content-based identifier for referencing lineage records.

- Provide a command-line interface for traversing and querying the lineage graph.

- Provide a way to compare two task runs to diagnose unexpected cache invalidation.

- Provide a way to verify that an output file has not been modified since it was produced.

## Non-goals

- Modifying the existing cache system -- lineage should be implemented as a supplemental layer to the cache and work directory.

- Committing to open standards such as CID, IPLD, OpenLineage -- the initial implementation should focus on supporting Nextflow-specific use cases, while leaving the door open to future integrations.

- Committing to a particular storage model -- the lineage system should expose an interface for *lineage stores* which can be implemented with any storage technology (filesystem, object storage, database, etc).

## Decision

Implement a lineage system in the Nextflow runtime that records every workflow run, task run, and output file as a standard *lineage record* with a unique *lineage ID* (LID). Implement a command-line interface for interacting with lineage records, including traversal and search.

## Core Capabilities

### Lineage records

Nextflow data lineage is modeled using the following types:

- `WorkflowRun`: represents a workflow execution (pipeline repository and revision, parameters, configuration)
- `TaskRun`: represents a task execution (script, container, inputs)
- `WorkflowOutput`: represents the output of a workflow run (based on the `output` block)
- `TaskOutput`: represents the output of a task run (based on the process `output:` section)
- `FileOutput`: represents an output file (file path, checksum, size, timestamps)

Every record has the same top-level structure:

- `version`: lineage schema version
- `kind`: lineage record type (enumerated above)
- `spec`: lineage record content

This approach is inspired by [Kubernetes API resources](https://kubernetes.io/docs/reference/using-api/api-concepts/). The initial lineage version is `v1beta1`, which will be replaced by `v1` when the lineage model is stable.

The [lineage schema](https://raw.githubusercontent.com/nextflow-io/schemas/main/lineage/v1beta1/schema.json) is available on GitHub.

### Lineage IDs (LIDs)

Every lineage record has a **lineage ID (LID)**. LIDs are unique and content-addressable.

- The `TaskRun` LID is equivalent to the existing Nextflow task hash. This is done primarily to make it easy for users to work with task lineage using the same task hashes that they already know (e.g. comparing two task runs using `nextflow lineage diff`).

- The `WorkflowRun` LID is derived from a hash of the lineage record itself.

- The LID of a `WorkflowOutput` / `TaskOutput` is the LID of the corresponding `WorkflowRun` / `TaskRun` suffixed with `#output`.

- The LID of a `FileOutput` is the LID of the `WorkflowRun` / `TaskRun` that produced the file, suffixed with the relative file path. For example, given a task with LID `05dfa13d` that produces the file `align/index.bam` in the task directory, the LID of the output file is `05dfa13d/align/index.bam`. For workflow outputs, the file path is relative to the workflow output directory (`-output-dir`).

Lineage records refer to each other by LID, forming an acyclic graph that can be traversed. For example:

- The `taskRun` field in `FileOutput` specifies the LID of the `TaskRun` that produced the output file
- The `input` field in `TaskRun` specifies the task inputs, including files which are specified as `FileOutput` LIDs (which in turn refer to the `TaskRun` that produced them)
- The `workflowRun` field in `TaskRun` specifies the LID of the workflow run that executed the task

LIDs can be referenced as URIs in Nextflow using the `lid://` URI scheme. This way, users can refer to files by their corresponding `FileOutput` LID instead of their absolute path:

```groovy
file('lid://05dfa13d/align/index.bam')
```

### Lineage store

A *lineage store* is a key-value store of lineage IDs to lineage records. The lineage system defines a simple interface for lineage stores, as well as a default implementation which stores records as JSON files in a directory tree (e.g. `.lineage/862df531.../.data.json`).

The lineage store is an extension point, so third-party plugins can integrate their own backends seamlessly into the core lineage system. For example, an alternative backend could store lineage records in a SQLite or Postgres database.

### Lineage generation

Lineage records are saved during the workflow run using a trace observer:

- When the workflow begins, save a `WorkflowRun` record
- When a task completes, save a `TaskRun` and `TaskOutput` record, as well as a `FileOutput` record for every output file
- When a file is published, save a `FileOutput` record linked to the workflow run
- When the workflow completes, save a `WorkflowOutput` record

Lineage generation must be explicitly enabled by setting `lineage.enabled = true` in the Nextflow configuration.

### Lineage CLI

The `nextflow lineage` command group provides several commands for interacting with data lineage:

- `list`: list workflow runs from the history log
- `view`: view a lineage record
- `find`: search for records matching an LID or labels
- `check`: check whether an output file has been modified since it was produced (based on the checksum in the `FileOutput` record)
- `diff`: compare two lineage records (useful for cache-invalidation diagnostics)
- `render`: render the lineage graph for a given record as an HTML diagram

The `view` and `find` commands always produce JSON output, making them easy to compose using JSON processing tools such as `jq`.

### Labels

Workflow output files can specify user-defined *labels* in the corresponding `FileOutput` record. Labels can be set via the `label` directive in the `output` block, and they can be used when searching for lineage records via `nextflow lineage find`.

Workflow output labels have limited granularity -- they can be set per-run or per-output, but not per-sample or per-file. They have been implemented as an experimental feature, but it is not yet clear what specific use cases will be satisfied by coarse-grained labels (experiment ID, QC vs alignment results, etc).

A workflow output can also produce an *index file*, which serializes the published channel as a structured file (e.g. CSV, JSON, YAML) containing published file paths and associated metadata. This metadata is effectively a fine-grained alternative to workflow output labels.

### Integration with Seqera Platform

While it remains to be seen how external systems such as Seqera Platform will surface data lineage to users, the lineage system has been designed with the following principles in mind:

- The lineage graph will be the source of truth for all data produced by Nextflow
- Searching lineage by metadata will be the primary way to find data
- An index is needed to facilitate real-time search (labels, output index files, etc)
- Datasets or *data products* can be implemented as saved lineage queries, allowing the data view to stay up-to-date as new data arrives

Searching a lineage graph using domain-specific metadata allows users to access their data directly, rather than going through runs or datasets produced by runs.

## Alternatives

### Content Identifier (CID)

[CID](https://github.com/multiformats/cid) is a self-describing content-addressed identifier used to enable content-addressable storage. While it was developed in the context of peer-to-peer networking protocols, it can be used on its own as a way to reference arbitrary data blobs based on their content.

Lineage IDs (LIDs) described above are inspired by the CID standard. LIDs are simpler and specialized for Nextflow, whereas CIDs are more generic because they are designed for interoperability between different data models. Using a Nextflow-specific identifier allows us to iterate rapidly based on our needs. However, we could adopt the CID standard as a future improvement.

### OpenLineage

[OpenLineage](https://openlineage.io/) is a specification for collecting and analyzing data lineage information. It is supported by several popular workflow systems such as Spark and Airflow.

One option for Nextflow was to simply adopt the OpenLineage standard for defining lineage records. This would allow Nextflow data lineage to be consumed by existing tools that support OpenLineage, such as [Marquez](https://marquezproject.ai/).

However, ultimately it made more sense to develop a Nextflow-specific lineage model, so that we can iterate and adapt the lineage model to our needs. OpenLineage could be supported in the future as an output format, i.e. adding the ability to export Nextflow lineage records to OpenLineage.

### Avoiding the use of `#output` in LIDs

Workflow outputs (`WorkflowOutput`) are recorded separately from workflow runs (`WorkflowRun`) because the identity of a run is determined only by the code and inputs, not the outputs. The same is true for task outputs (`TaskOutput`) and runs (`TaskRun`).

This approach presents a dilemma. The output must be recorded separately, because an upstream record cannot depend on a downstream record (the lineage graph must be *acyclic*). However, this means that a workflow output cannot be retrieved directly because it is downstream of the workflow run.

The current solution is to define a pseudo-LID for workflow outputs as the workflow run LID suffixed with `#output`. This way, the workflow output can be accessed directly without compromising the identity of the `WorkflowRun` record or the acyclic nature of the lineage graph.

An alternative approach is to refactor the workflow records as follows:

- `WorkflowLaunch`: created on workflow launch (replaces `WorkflowRun`)
- `WorkflowRun`: created on workflow completion, contains a reference to `WorkflowLaunch` record and the contents of `WorkflowOutput`

This approach would eliminate the need for the `#output` suffix. Instead, the `WorkflowRun` record would be saved when the run completes rather than when it starts. It also matches the model of runs vs launch configurations in Seqera Platform. The same approach could be applied to task runs and outputs.

The `#output` suffix is part of a broader strategy to align the lineage system with the existing runtime:

- Use the existing task hash as the `TaskRun` LID
- Use `#output` to reference workflow / task outputs from the run LID
- Use the relative file path to reference `FileOutput` records from the run LID
- Try to use the same conventions for workflow runs and task runs

A "pure LID" approach that is closer to the CID model remains an option for future improvement, if it is deemed worthwhile. However, the current approach is practical and allows us to introduce lineage as a supplemental layer without modifying the core runtime.

See also: [#6011](https://github.com/nextflow-io/nextflow/pull/6011)

### Unifying the lineage, cache, and work directory

Nextflow data lineage is implemented as a supplemental layer on the existing data model (metadata cache, work directory). Lineage records refer to output files by their work directory path. As a result, lineage records are tied to a specific storage location -- if the work directory were moved, the corresponding lineage records would be invalidated.

One alternative is to unify the lineage, metadata cache, and work directory into a single data model. For example, task outputs could be stored directly in the lineage store as content-addressable blobs instead of files in the work directory. This would be more coherent overall, and as a side benefit it would be easier to move the data around.

However, the current approach is useful for experimenting with lineage without affecting the existing cache and work directory, which are core components of the Nextflow runtime. It should be investigated in the future when the lineage system is more mature.

## Links

- [nf-core/bytesize: Content Addressable Data Storage](https://www.youtube.com/watch?v=7wmaWU6-pYM)
- [nf-blocks](https://github.com/robsyme/nf-blocks) plugin
- Related: [Workflow outputs ADR](./20251020-workflow-outputs.md)
