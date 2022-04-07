/*
 * Copyright (c) 2020-2021. Seqera Labs, S.L.
 *
 * All Rights reserved
 *
 */

package nextflow.cloud.google.batch.model

import groovy.transform.CompileStatic
/**
 * Model a Batch Task Group
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class TaskGroup {
    TaskSpec taskSpec
    Map<String,String> labels

    TaskGroup withTaskSpec(TaskSpec spec) {
        this.taskSpec = spec
        return this
    }

    TaskGroup withLabels(Map<String,String> labels) {
        this.labels = labels
        return this
    }
}
