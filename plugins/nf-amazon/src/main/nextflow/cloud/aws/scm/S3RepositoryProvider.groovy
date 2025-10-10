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
import org.eclipse.jgit.lib.Constants
import org.eclipse.jgit.lib.Ref
import org.eclipse.jgit.transport.CredentialsProvider

import java.nio.file.Files
import java.nio.file.Path


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

    /**
     * {@inheritDoc}
     *
     * Note: S3 git-remote stores repositories as Git bundles in S3 (one bundle per branch).
     * Reading a single file requires downloading and unpacking the entire bundle for that branch.
     * When no revision is specified, we determine the default branch from the remote HEAD
     * to avoid downloading unnecessary branches.
     */
    @Override
    byte[] readBytes( String path ) {
        log.debug("Reading $path from S3 git-remote")
        Path tmpDir = null
        try {
            tmpDir = Files.createTempDirectory("s3-git-remote-")

            // Determine which branch to clone
            def branchToClone = revision
            if (!branchToClone) {
                // No revision specified - fetch only the default branch
                // This avoids downloading unnecessary branch bundles
                branchToClone = getDefaultBranch()
                log.debug("No revision specified, using default branch: $branchToClone")
            }

            final command = Git.cloneRepository()
                .setURI(getEndpointUrl())
                .setDirectory(tmpDir.toFile())
                .setCredentialsProvider(getGitCredentials())
                .setCloneAllBranches(false)  // Only clone the specified branch
                .setBranch(branchToClone)

            command.call()
            final file = tmpDir.resolve(path)
            return Files.exists(file) ? Files.readAllBytes(file) : null
        }
        catch (Exception e) {
            log.debug("Unable to retrieve file: $path from repo: $project", e)
            return null
        }
        finally {
            if (tmpDir != null && Files.exists(tmpDir)) {
                tmpDir.toFile().deleteDir()
            }
        }
    }

    /**
     * Get the default branch from the S3 git-remote repository by querying remote refs.
     * Uses Git's lsRemote to fetch the HEAD symbolic ref, which points to the default branch.
     *
     * @return The default branch name
     */
    @Memoized
    String getDefaultBranch() {
        // Fetch remote refs using Git's lsRemote
        final refs = fetchRefs()
        if (!refs){
            throw new Exception("No remote references found")
        }
        // Find the HEAD symbolic ref
        final headRef = refs.find { it.name == Constants.HEAD }

        if (!headRef )
            throw new Exception("No remote HEAD ref found ")

        if( !headRef.isSymbolic() )
            throw new Exception("Incorrect HEAD ref. Not a symbolic ref.")

        final target = headRef.target.name
        if( target.startsWith('refs/heads/') ) {
            return target.substring('refs/heads/'.length())
        }
        return target
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
