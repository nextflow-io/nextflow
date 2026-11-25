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

package nextflow.module

import java.util.regex.Pattern

import groovy.transform.CompileStatic
import groovy.transform.Immutable
import io.seqera.npr.api.schema.v1.Pipeline
import io.seqera.npr.api.schema.v1.PipelineRelease
import io.seqera.npr.api.schema.v1.PipelineReleaseStatus
import io.seqera.npr.client.RegistryClient
import io.seqera.npr.client.RegistryException
import nextflow.exception.AbortOperationException

/**
 * Release of a registry pipeline: where to fetch it from and the registry it comes from
 */
@Immutable
@CompileStatic
class PipelineSource {

    private static final Pattern NAME_PATTERN = ~/^[a-z0-9][a-z0-9._-]*\/[a-z][a-z0-9._-]*$/

    String name
    String version
    String workflowPath
    String repositoryUrl
    String revision
    String commitId
    String registryUrl

    static boolean isValidName(String name) {
        return name && NAME_PATTERN.matcher(name).matches()
    }

    /**
     * Release of a registry pipeline, the latest when no version is given; fails when it is missing or revoked
     */
    static PipelineSource resolve(RegistryClient client, String name, String version) {
        if( !isValidName(name) )
            throw new AbortOperationException("Invalid pipeline name '$name' -- expected <scope>/<name>, e.g. nf-core/rnaseq")
        final Pipeline pipeline
        final PipelineRelease release
        try {
            pipeline = client.getPipeline(name)
            release = version ? client.getPipelineRelease(name, version) : pipeline.latest
        }
        catch( RegistryException e ) {
            if( e.statusCode == 404 )
                throw new AbortOperationException("Pipeline ${name}${version ? '@' + version : ''} not found in the registry", e)
            throw e
        }
        if( release == null )
            throw new AbortOperationException("Pipeline ${name} has no published release")
        if( pipeline.revoked || release.status == PipelineReleaseStatus.REVOKED )
            throw new AbortOperationException("Pipeline ${name}@${release.version} is revoked")
        return new PipelineSource(
            name: name,
            version: release.version,
            workflowPath: pipeline.workflowPath,
            repositoryUrl: pipeline.repositoryUrl,
            revision: release.gitRevision,
            commitId: release.commitId,
            registryUrl: client.primaryUrl())
    }
}
