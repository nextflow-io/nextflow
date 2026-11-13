/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.util

import groovy.transform.CompileStatic
import nextflow.script.WorkflowMetadata

import java.nio.file.Path

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class PathNormalizer {

    private URL repository

    private String commitId

    private String projectDir

    private String workDir

    PathNormalizer(WorkflowMetadata metadata) {
        repository = metadata.repository ? new URL(metadata.repository) : null
        commitId = metadata.commitId
        projectDir = metadata.projectDir.normalize().toUriString()
        workDir = metadata.workDir.normalize().toUriString()
    }

    /**
     * Normalize paths against the original remote URL, or
     * work directory, where appropriate.
     *
     * @param path
     */
    String normalizePath(Path path) {
        normalizePath(path.toUriString())
    }

    String normalizePath(String path) {
        if(!path)
            return null
        // replace work directory with relative path
        if( path.startsWith(workDir) )
            return path.replace(workDir, 'work')

        // replace project directory with source URL (if applicable)
        if( repository && path.startsWith(projectDir) )
            return getProjectSourceUrl(path)

        // encode local absolute paths as file URLs
        if( path.startsWith('/') )
            return 'file://' + path

        return path
    }

    /**
     * Get the source URL for a project asset.
     *
     * @param path
     */
    private String getProjectSourceUrl(String path) {
        switch( repository.host ) {
        case 'bitbucket.org':
            return path.replace(projectDir, "${repository}/src/${commitId}")
        case 'github.com':
            return path.replace(projectDir, "${repository}/tree/${commitId}")
        case 'gitlab.com':
            return path.replace(projectDir, "${repository}/-/tree/${commitId}")
        default:
            return path
        }
    }

}
