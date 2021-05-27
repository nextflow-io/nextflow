package io.nextflow.gradle


import groovy.transform.CompileStatic
import io.nextflow.gradle.model.S3Extension
import org.gradle.api.Plugin
import org.gradle.api.Project
/**
 * Implements Nextflow build plugin
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */

@CompileStatic
class NextflowBuildPlugin implements Plugin<Project> {

    void apply(Project target) {
        target.extensions.create('s3', S3Extension)
    }

}

