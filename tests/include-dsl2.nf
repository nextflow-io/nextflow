#!/bin/bash nextflow

// Original way of module inclusion with an alias in DSL2
// Note, "test1" can be a GString here, but test1_alias cannot.
include {test1 as test1_alias;} from './subworkflow-dsl2.nf'

// Enhanced functionality: ability to use a GString when setting an included
// module's alias:
workflow="test2"
workflow_alias="test2_alias"
include {"${workflow}" as "${workflow_alias}";} from './subworkflow-dsl2.nf'

workflow {
    test1_alias()
    "${workflow_alias}"()
}
