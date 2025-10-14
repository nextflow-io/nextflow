/*
 * Copyright 2024-2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.config

import groovy.transform.CompileStatic
import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

@ScopeName("workflow")
@Description("""
    The `workflow` scope provides workflow execution options.
""")
@CompileStatic
class WorkflowConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        When `true`, the pipeline will exit with a non-zero exit code if any failed tasks are ignored using the `ignore` error strategy (default: `false`).
    """)
    final boolean failOnIgnore

    @ConfigOption
    @Description("""
        Specify a closure that will be invoked at the end of a workflow run (including failed runs).
    """)
    final Closure onComplete

    @ConfigOption
    @Description("""
        Specify a closure that will be invoked if a workflow run is terminated.
    """)
    final Closure onError

    @Description("""
        The `workflow.output` scope provides options for publishing workflow outputs.
    
        [Read more](https://nextflow.io/docs/latest/reference/config.html#workflow)
    """)
    final WorkflowOutputConfig output

    /* required by extension point -- do not remove */
    WorkflowConfig() {}

    WorkflowConfig(Map opts) {
        failOnIgnore = opts.failOnIgnore as boolean
        onComplete = opts.onComplete as Closure
        onError = opts.onError as Closure
        output = new WorkflowOutputConfig(opts.output as Map ?: Collections.emptyMap())
    }

}

@CompileStatic
class WorkflowOutputConfig implements ConfigScope {

    @ConfigOption
    @Description("""
        *Currently only supported for S3.*

        Specify the media type, also known as [MIME type](https://developer.mozilla.org/en-US/docs/Web/HTTP/MIME_types), of published files (default: `false`). Can be a string (e.g. `'text/html'`), or `true` to infer the content type from the file extension.
    """)
    final Object contentType

    @ConfigOption
    @Description("""
        *Currently only supported for local and shared filesystems.*

        Copy file attributes (such as the last modified timestamp) to the published file (default: `false`).
    """)
    final boolean copyAttributes

    @ConfigOption
    @Description("""
        Enable or disable publishing (default: `true`).
    """)
    final boolean enabled

    @ConfigOption
    @Description("""
        When `true`, the workflow will not fail if a file can't be published for some reason (default: `false`).
    """)
    final boolean ignoreErrors

    @ConfigOption
    @Description("""
        The file publishing method (default: `'symlink'`).
    """)
    final String mode

    @ConfigOption
    @Description("""
        When `true` any existing file in the specified folder will be overwritten (default: `'standard'`).
    """)
    final Object overwrite

    @ConfigOption
    @Description("""
        *Currently only supported for S3.*

        Specify the storage class for published files.
    """)
    final String storageClass

    @ConfigOption
    @Description("""
        *Currently only supported for S3.*

        Specify arbitrary tags for published files.
    """)
    final Map tags

    /* required by extension point -- do not remove */
    WorkflowOutputConfig() {}

    WorkflowOutputConfig(Map opts) {
        contentType = opts.contentType
        copyAttributes = opts.copyAttributes as boolean
        enabled = opts.enabled as boolean
        ignoreErrors = opts.ignoreErrors as boolean
        mode = opts.mode
        overwrite = opts.overwrite
        storageClass = opts.storageClass
        tags = opts.tags as Map
    }

}
