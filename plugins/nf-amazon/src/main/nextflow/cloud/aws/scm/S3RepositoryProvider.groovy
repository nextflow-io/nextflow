/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.cloud.aws.scm

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.cloud.aws.scm.jgit.S3GitCredentialsProvider
import nextflow.exception.AbortOperationException
import nextflow.scm.ProviderConfig
import nextflow.scm.RepositoryProvider
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.api.errors.TransportException
import org.eclipse.jgit.transport.CredentialsProvider

import java.nio.file.Files


/**
 * Implements a repository provider for git-remote-s3 repositories.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class S3RepositoryProvider extends RepositoryProvider {

    S3RepositoryProvider(String project, ProviderConfig config) {
        assert config instanceof S3ProviderConfig
        log.debug("Creating S3 repository provider for $project")
        this.project = project
        this.config = config
    }
    /** {@inheritDoc} **/
    @Memoized
    @Override
    CredentialsProvider getGitCredentials() {
        final providerConfig = this.config as S3ProviderConfig
        final credentials = new S3GitCredentialsProvider()
        if( providerConfig.region )
            credentials.setRegion(providerConfig.region)
        if( providerConfig.awsCredentialsProvider )
            credentials.setAwsCredentialsProvider(providerConfig.awsCredentialsProvider)
        return credentials
    }

    /** {@inheritDoc} **/
    // called by AssetManager
    // used to set credentials for a clone, pull, fetch, operation
    @Override
    boolean hasCredentials() {
        // set to true
        // uses AWS Credentials instead of username : password
        // see getGitCredentials()
        return true
    }

    /** {@inheritDoc} **/
    @Override
    String getName() { return project }

    /** {@inheritDoc} **/
    @Override
    String getEndpointUrl() {
        return "s3://$project"
    }

    /** {@inheritDoc} **/
    // not used, but the abstract method needs to be overridden
    @Override
    String getContentUrl( String path ) {
        throw new UnsupportedOperationException()
    }

    /** {@inheritDoc} **/
    // called by AssetManager
    @Override
    String getCloneUrl() { getEndpointUrl() }

    /** {@inheritDoc} **/
    // called by AssetManager
    @Override
    String getRepositoryUrl() { getEndpointUrl() }

    /** {@inheritDoc} **/
    // called by AssetManager
    // called by RepositoryProvider.readText()
    @Override
    byte[] readBytes( String path ) {
        log.debug("Reading $path")
        //Not possible to get a single file requires to clone the branch and get the file
        final tmpDir = Files.createTempDirectory("s3-git-remote")
        final command = Git.cloneRepository()
            .setURI(getEndpointUrl())
            .setDirectory(tmpDir.toFile())
            .setCredentialsProvider(getGitCredentials())
        if( revision )
            command.setBranch(revision)
        try {
            command.call()
            final file = tmpDir.resolve(path)
            return file.getBytes()
        }
        catch (Exception e) {
            log.debug(" unable to retrieve file: $path from repo: $project", e)
            return null
        }
        finally{
            tmpDir.deleteDir()
        }
    }

    @Override
    List<RepositoryEntry> listDirectory(String path, int depth) {
        throw new UnsupportedOperationException("S3-git-remote does not support 'listDirectory' operation")
    }

/** {@inheritDoc} **/
    // called by AssetManager
    @Override
    void validateRepo() {
        // Nothing to check
    }

    private String errMsg(Exception e) {
        def msg = "Unable to access Git repository"
        if( e.message )
            msg + " - ${e.message}"
        else
            msg += ": " + getCloneUrl()
        return msg
    }

    @Override
    List<BranchInfo> getBranches() {
        try {
            return super.getBranches()
        }
        catch ( TransportException e) {
            throw new AbortOperationException(errMsg(e), e)
        }
    }

    @Override
    List<TagInfo> getTags() {
        try {
            return super.getTags()
        }
        catch (TransportException e) {
            throw new AbortOperationException(errMsg(e), e)
        }
    }


}
