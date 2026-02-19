(data-lineage-page)=

# Getting started with data lineage

Data lineage in Nextflow provides comprehensive tracking of workflow runs, task executions, and output files. This feature helps you verify the integrity and reproducibility of your pipeline results by maintaining a complete history of computations and intermediate data.

:::{warning}
Data lineage is an experimental feature added in Nextflow 25.04. The functionality may change in future releases.
:::

## Overview

Data lineage tracks the complete history of your Nextflow runs, including:

- Workflow runs and their configurations
- Task executions and their inputs/outputs
- File outputs and their provenance

Each lineage record has a unique identifier called a *lineage ID* (LID) that you can use to access and query the data.

:::{note}
The data model for every lineage record is defined in the Nextflow [source code](https://github.com/nextflow-io/nextflow/tree/master/modules/nf-lineage/src/main/nextflow/lineage/model).
:::

## Enable data lineage

To enable data lineage tracking, add the following to your Nextflow configuration:

```groovy
lineage.enabled = true
```

By default, lineage data is stored in the `.lineage` directory in your current working directory. You can customize this location:

```groovy
lineage.store.location = '<PATH_TO_STORAGE>'
```

:::{tip}
For global configuration, add these settings to `$HOME/.nextflow/config`.
:::

The location can be local or remote (e.g., an s3 storage bucket). See the {ref}`config-lineage` configuration scope for details.

## Generate lineage metadata

Run a Nextflow pipeline to generate some lineage metadata. For example:

```console
$ nextflow run rnaseq-nf -profile conda
```

Nextflow will automatically record the workflow run, task executions, and output files in the lineage store.

## Explore lineage

Now that you have generated some lineage metadata, you can explore it from the command line using the {ref}`cli-lineage` command.

First, use the `list` subcommand to list the workflow runs in the lineage store:

```console
$ nextflow lineage list
TIMESTAMP               RUN NAME                SESSION ID                              LINEAGE ID                            
2025-05-09 13:28:30 CDT peaceful_blackwell      065bdc6b-89b4-42ee-92c1-2a5af37f2c50    lid://16b31030474f2e96c55f4940bca3ab64
```

The *lineage ID* (LID) is the unique identifier for the workflow run and the entrypoint for exploring the lineage.

Use the `view` subcommand to view the lineage record for the workflow run:

```console
$ nextflow lineage view lid://16b31030474f2e96c55f4940bca3ab64
{
  "version": "lineage/v1beta1",
  "type": "WorkflowRun",
  "workflow": {
    "scriptFiles": [
      ...
    ],
    "repository": "https://github.com/nextflow-io/rnaseq-nf",
    "commitId": "86165b8c81d43a1f57363964431395152e353e56"
  },
  "sessionId": "065bdc6b-89b4-42ee-92c1-2a5af37f2c50",
  "name": "peaceful_blackwell",
  "params": [
    ...
  ],
  "config": {
    ...
  }
}
```

Every workflow run is represented in the lineage store as a `WorkflowRun` record. It includes information such as the pipeline repository, revision, run name, parameters, and resolved config.

The output files of a workflow run can be accessed as `lid://<WORKFLOW_RUN_HASH>/<PATH>`, where `<PATH>` is the file path relative to the workflow output directory.

:::{note}
Files must be published to the workflow output directory as defined by the `outputDir` config option (or `-output-dir` command line option) in order to be recorded as workflow outputs in the lineage store.
:::

List the output directory to see the available files:

```console
$ find results
results
results/fastqc_ggal_gut_logs
results/fastqc_ggal_gut_logs/ggal_gut_1_fastqc.html
results/fastqc_ggal_gut_logs/ggal_gut_1_fastqc.zip
results/fastqc_ggal_gut_logs/ggal_gut_2_fastqc.html
results/fastqc_ggal_gut_logs/ggal_gut_2_fastqc.zip
results/multiqc_report.html
```

Now, use the workflow LID and relative path to view the lineage record for an output file:

```console
$ nextflow lineage view lid://16b31030474f2e96c55f4940bca3ab64/multiqc_report.html
{
  "version": "lineage/v1beta1",
  "type": "FileOutput",
  "path": "/results/multiqc_report.html",
  "checksum": {
    "value": "03fd5ed150c7862e1fad5efd4f574a47",
    "algorithm": "nextflow",
    "mode": "standard"
  },
  "source": "lid://862df53160e07cd823c0c3960545e747/multiqc_report.html",
  "workflowRun": "lid://16b31030474f2e96c55f4940bca3ab64",
  "taskRun": null,
  "size": 5079806,
  "createdAt": "2025-05-09T13:27:34.576590545-05:00",
  "modifiedAt": "2025-05-09T13:27:34.586590551-05:00",
  "labels": null
}
```

Every output file is represented in the lineage store as a `FileOutput` record. It includes basic file information, such as the real path, checksum, file size and created/modified timestamps, as well as lineage information, such as the workflow run and task run that produced it.

As this record is a workflow output, it is not linked directly to a task run. Instead, it is linked to the original task output.

Any LID in a lineage record can be viewed, allowing you to traverse the lineage metadata interactively. Use the value of `source` to view the original task output:

```console
$ nextflow lineage view lid://862df53160e07cd823c0c3960545e747/multiqc_report.html
{
  "version": "lineage/v1beta1",
  "type": "FileOutput",
  "path": "/work/86/2df53160e07cd823c0c3960545e747/multiqc_report.html",
  "checksum": {
    "value": "b14f5171a48ce5c22ea27d7b8e57b6c4",
    "algorithm": "nextflow",
    "mode": "standard"
  },
  "source": "lid://862df53160e07cd823c0c3960545e747",
  "workflowRun": "lid://16b31030474f2e96c55f4940bca3ab64",
  "taskRun": "lid://862df53160e07cd823c0c3960545e747",
  "size": 5079806,
  "createdAt": "2025-05-09T13:27:34.236590379-05:00",
  "modifiedAt": "2025-05-09T13:27:34.246590383-05:00",
  "labels": null
}
```

This record is the task output for the same file -- it has a value for `taskRun` which is the same as its `source`.

View the lineage record for the task that produced this file:

```console
$ nextflow lineage view lid://862df53160e07cd823c0c3960545e747
{
  "version": "lineage/v1beta1",
  "type": "TaskRun",
  "sessionId": "065bdc6b-89b4-42ee-92c1-2a5af37f2c50",
  "name": "MULTIQC",
  "codeChecksum": {
    "value": "edf2e9f84cd3a18ee9259012b660f2dd",
    "algorithm": "nextflow",
    "mode": "standard"
  },
  "script": "\n    cp multiqc/* .\n    echo \"custom_logo: $PWD/nextflow_logo.png\" \u003e\u003e multiqc_config.yaml\n    multiqc -n multiqc_report.html .\n    ",
  "input": [
    {
      "type": "path",
      "name": "*",
      "value": [
        "lid://eff8846883b46c5a76f11e7e4480a6c8/ggal_gut",
        "lid://2d8bd92c69f732605bc99941e60d5319/fastqc_ggal_gut_logs"
      ]
    },
    {
      "type": "path",
      "name": "config",
      "value": [
        {
          "path": "https://github.com/nextflow-io/rnaseq-nf/tree/86165b8c81d43a1f57363964431395152e353e56/multiqc",
          "checksum": {
            "value": "2aac500cdfb292e961e678433e7dc3d8",
            "algorithm": "nextflow",
            "mode": "standard"
          }
        }
      ]
    }
  ],
  "container": null,
  "conda": "file:///conda/env-4a436c230263dfdbbf4dddd0623505d1",
  "spack": null,
  "architecture": null,
  "globalVars": {},
  "binEntries": [],
  "workflowRun": "lid://16b31030474f2e96c55f4940bca3ab64"
}
```

Every task run is represented in the lineage store as a `TaskRun`, which includes information such as the name, script, inputs, and software dependencies. From here, you can continue traversing through the file inputs to view upstream tasks.

Finally, use the `render` subcommand to render the entire lineage of the MULTIQC report as an HTML report:

```console
$ nextflow lineage render lid://16b31030474f2e96c55f4940bca3ab64/multiqc_report.html
Rendered lineage graph for lid://16b31030474f2e96c55f4940bca3ab64/multiqc_report.html to lineage.html
```

Open the HTML report in a web browser to view the lineage graph.

## Query lineage records

To find a lineage record, you normally have to know the LID of the record or a downstream record (such as a workflow run) from which you can traverse to the desired record. However, you can also query the entire lineage store by fields to quickly find relevant records and aggregate records from different runs.

Use the `find` subcommand to find all tasks executed by a workflow run:

```console
$ nextflow lineage find type=TaskRun workflowRun=lid://16b31030474f2e96c55f4940bca3ab64
[
  "lid://2d8bd92c69f732605bc99941e60d5319",
  "lid://eff8846883b46c5a76f11e7e4480a6c8",
  "lid://862df53160e07cd823c0c3960545e747",
  "lid://6d3bff36bf2c3c14c2d383384621e8ca"
]
```

You can use any field defined in the [lineage data model](https://github.com/nextflow-io/nextflow/tree/master/modules/nf-lineage/src/main/nextflow/lineage/model).

:::{tip}
Since the `find` and `view` subcommands always output JSON, you can use JSON processing tools such as [jq](https://jqlang.org/) to further query and transform results.
:::

## Compare task runs

Task run LIDs are based on the standard {ref}`task hash <cache-resume-task-hash>`, which makes it easy to compare two task runs in the lineage metadata. For example, if a task is unexpectedly re-executed during a resumed run, as long as lineage is enabled for both the initial and resumed runs, the two tasks can be compared without any additional runs.

This section builds on the above [`rnaseq-nf` example](#generate-lineage-metadata) to demonstrate how to compare two task runs in the event of a cache invalidation.

First, modify the pipeline in a way that invalidates the cache for the `MULTIQC` process. For example, modify the process script.

Resume the pipeline. It will re-execute the `MULTIQC` process:

```console
$ nextflow run rnaseq-nf -profile conda -resume

...

[6d/3bff36] process > RNASEQ:INDEX (ggal_1_48850000_49020000) [100%] 1 of 1, cached: 1 ✔
[2d/8bd92c] process > RNASEQ:FASTQC (FASTQC on ggal_gut)      [100%] 1 of 1, cached: 1 ✔
[ef/f88468] process > RNASEQ:QUANT (ggal_gut)                 [100%] 1 of 1, cached: 1 ✔
[94/33dda7] process > MULTIQC                                 [100%] 1 of 1 ✔
```

Retrieve the hash of the `MULTIQC` run from the log file or work directory. Compare it to the task hash of the initial run:

```console
$ nextflow lineage diff lid://862df53160e07cd823c0c3960545e747 lid://9433dda73f2193491f9a26e3e23cd8a1
diff --git 862df53160e07cd823c0c3960545e747 9433dda73f2193491f9a26e3e23cd8a1
--- 862df53160e07cd823c0c3960545e747
+++ 9433dda73f2193491f9a26e3e23cd8a1
@@ -3,11 +3,11 @@
   "sessionId": "065bdc6b-89b4-42ee-92c1-2a5af37f2c50",
   "name": "MULTIQC",
   "codeChecksum": {
-    "value": "edf2e9f84cd3a18ee9259012b660f2dd",
+    "value": "9615a8da3a3f9e935cfc8e4042cdf5e0",
     "algorithm": "nextflow",
     "mode": "standard"
   },
-  "script": "\n    cp multiqc/* .\n    echo \"custom_logo: $PWD/nextflow_logo.png\" \u003e\u003e multiqc_config.yaml\n    multiqc -n multiqc_report.html .\n    ",
+  "script": "\n    cp multiqc/* . # hello!\n    echo \"custom_logo: $PWD/nextflow_logo.png\" \u003e\u003e multiqc_config.yaml\n    multiqc -n multiqc_report.html .\n    ",
   "input": [
     {
       "type": "path",
@@ -38,5 +38,5 @@
   "architecture": null,
   "globalVars": {},
   "binEntries": [],
-  "workflowRun": "lid://16b31030474f2e96c55f4940bca3ab64"
+  "workflowRun": "lid://65044872aad36f97e42336b9ba0dee57"
 }
```

Note the difference between the task scripts, highlighting the change that caused the task to be re-executed.

## Use lineage with workflow outputs

Workflow outputs declared in the `output` block are also recorded in the lineage store. The output of a workflow run is accessible as `lid://<WORKFLOW_RUN_HASH>#output`.

For example, run the `rnaseq-nf` pipeline using the `preview-25-04` branch, which uses the `output` block to publish outputs:

```console
$ nextflow -r preview-25-04 -profile conda
```

View the workflow output in the lineage metadata:

```console
$ nextflow lineage view lid://9410d13abeec617640b5fe9735ba12fc#output
[
  {
    "type": "Collection",
    "name": "samples",
    "value": "lid://9410d13abeec617640b5fe9735ba12fc/samples.json"
  },
  {
    "type": "Path",
    "name": "summary",
    "value": "lid://9410d13abeec617640b5fe9735ba12fc/multiqc_report.html"
  }
]
```

This view can be used to traverse output files directly instead of inferring LIDs from the workflow output directory.

See {ref}`workflow-output-def` for more information about the `output` block.

## Use lineage in a Nextflow script

Since lineage IDs are valid URIs, output files in the lineage store can be accessed by their LID in a Nextflow script, like any other path. The LID path returns the *real* path as defined by the `path` field in the `FileOutput` record.

The following script uses the `samples.json` from the previous example as an input samplesheet:

```nextflow
channel.fromPath('lid://9410d13abeec617640b5fe9735ba12fc/samples.json')
    .splitJson()
    .view()
```

It should produce the following output:

```console
[id:gut, quant:/results/gut/quant, fastqc:/results/gut/fastqc]
```

The `fromLineage` channel factory can also be used to query lineage records in a similar manner as the `find` subcommand. See {ref}`channel-from-lineage` for details.
