
Introduction
============

The purpose of this document is to describe an idea for decomposing Nextflow workflows into subworkflows and then present a preliminary design for achieving this.

The overarching goal is to allow a Nextflow workflow to be broken into smaller "modules" for organizational purposes and to allow these "modules" to re-used within other workflows.  All of the arguments that hold for software modularity in general, hold here as well.  This document describes the idea where a Nextflow _workflow_ itself becomes a module, so there is really no concept of modules, but rather just workflows calling workflows with the same execution space.

Design
======

The central idea is to add two new collections to the Nextflow _workflow_ object created for every Nextflow script, a collection of input channels and a collection of output channels.  These input and output collections of channels would serve as the _interface_ to the workflow. If all input and output channels are defined in terms of `workflow.input` and `workflow.output`, then the workflow could be used as a module.

Depending on how a workflow script is invoked, the `workflow.input` field may be full or empty. Invoking a workflow in standalone mode would result in an empty `workflow.input` and the assumption is that all channels would need to be populated as they are now.  Invoking a workflow as a `subworkflow` would populate `workflow.input` and would allow channels defined in another workflow script to be passed in.

Subworkflows will be consumed by other workflows by extending `process` with a new component called _subworkflow_, which will reference the subworkflow.  The `input:` and `output:` channels of the process will map to the `workflow.input` and `workflow.output` channels defined in the subworkflow.


Example
======

```groovy

Channel
    .fromPath("wherever")
    .set{ a_channel }

process normal {
    input:
    file(a) from a_channel //

    output:
    stdout into x_channel

    //...
}


process coolthing {

    // Input channel names must line up with subworkflow's workflow.input.
    input:
    input_1 from x_channel

    // Output channels must line up with subworkflow's workflow.output.
    output:
    output_x into x_results

    // Behind the scenes, nextflow loads the workflow script identified in this
    // subworkflow section and maps the channels identified above into the workflow.input
    // and workflow.output of the subworkflow script.
    subworkflow:
      path: 'sub/coolthing_module.nf'

//// or...

    subworkflow:
      git: 'http://coolthing/nextflow/repo'
      tag: 'v2.3.4'
}

process usual {

    input:
    val(x) from x_results

    output:
    stdout into whatever

    //...

}
```

Here is one possible syntax for a _subworkflow_ file:

```groovy

// Allow the script to be used standalone, where the input channel
// is defined normally if workflow.input is empty.
if (workflow.input.isEmpty()) {
    input_1 = Channel.fromPath(params.whatever)

// Otherwise, if workflows.input is populated, lookup the desired
// channels and use them.
} else {
    input_1 = workflow.input["input_1"]
}

process coolthing_a {
    input:
    val(x) from input_1

    output:
    stdout into output_1

    //...
}

process coolthing_b {
    input:
    val(xx) from output_1

    // Output channels are always put into the workflow.output collection.
    output:
    stdout into workflow.output["output_x"]

    //...
}
```

