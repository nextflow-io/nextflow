/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.cloud.aws.codecommit

import software.amazon.awssdk.auth.credentials.AwsBasicCredentials
import software.amazon.awssdk.auth.credentials.DefaultCredentialsProvider
import software.amazon.awssdk.auth.credentials.StaticCredentialsProvider
import software.amazon.awssdk.core.exception.SdkException
import software.amazon.awssdk.regions.Region
import software.amazon.awssdk.services.codecommit.CodeCommitClient
import software.amazon.awssdk.services.codecommit.model.CodeCommitException
import software.amazon.awssdk.services.codecommit.model.GetFileRequest
import software.amazon.awssdk.services.codecommit.model.GetFolderRequest
import software.amazon.awssdk.services.codecommit.model.GetRepositoryRequest
import software.amazon.awssdk.services.codecommit.model.RepositoryMetadata
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.exception.MissingCredentialsException
import nextflow.scm.ProviderConfig
import nextflow.scm.RepositoryProvider
import nextflow.scm.RepositoryProvider.RepositoryEntry
import nextflow.util.StringUtils
import org.eclipse.jgit.api.errors.TransportException
import org.eclipse.jgit.transport.CredentialsProvider
/**
 * Implements a repository provider for AWS CodeCommit
 *
 * @author W. Lee Pang <wleepang@gmail.com>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AwsCodeCommitRepositoryProvider extends RepositoryProvider {

    AwsCodeCommitRepositoryProvider(String project, ProviderConfig config) {
        assert config instanceof AwsCodeCommitProviderConfig
        this.project = project  // expect: "codecommit-<region>/<repository>"
        this.config = config
        this.region = config.region
        this.repositoryName = project.tokenize('/')[-1]
        this.client = createClient(config)
    }

    private String region
    private CodeCommitClient client
    private String repositoryName


    protected CodeCommitClient createClient(AwsCodeCommitProviderConfig config) {
        final builder = CodeCommitClient.builder()
                .region(Region.of(region))
        if( config.user && config.password ) {
            final creds = AwsBasicCredentials.create(config.user, config.password)
            log.debug "AWS CodeCommit using username=$config.user; password=${StringUtils.redact(config.password)}"
            builder.credentialsProvider( StaticCredentialsProvider.create(creds) )
        }
        else {
            log.debug "AWS CodeCommit using default credentials chain"
            builder.credentialsProvider( DefaultCredentialsProvider.builder().build() )
        }
        builder.build()
    }

    /** {@inheritDoc} **/
    @Memoized
    @Override
    CredentialsProvider getGitCredentials() {
        return new AwsCodeCommitCredentialProvider(username: user, password: password)
    }

    private RepositoryMetadata getRepositoryMetadata() {
        final request = GetRepositoryRequest.builder()
                .repositoryName(repositoryName)
                .build()

        return client
                .getRepository(request)
                .repositoryMetadata()
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
    String getName() { "CodeCommit" }

    /** {@inheritDoc} **/
    @Override
    String getEndpointUrl() {
        "https://git-codecommit.${region}.amazonaws.com/v1/repos/${repositoryName}"
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

        final builder = GetFileRequest.builder()
            .repositoryName(repositoryName)
            .filePath(path)
        if( revision )
            builder.commitSpecifier(revision)

        try {
            return client
                    .getFile( builder.build() )
                    .fileContent()?.asByteArray()
        }
        catch (Exception e) {
            checkMissingCredsException(e)
            log.debug "AWS CodeCommit unable to retrieve file: $path from repo: $repositoryName"
            return null
        }
    }

    /** {@inheritDoc} **/
    @Override
    List<RepositoryEntry> listDirectory(String path, int depth) {
        try {
            // AWS CodeCommit doesn't have a dedicated directory listing API like GitHub
            // We would need to use GetFolder API, but it has limitations
            def request = GetFolderRequest.builder()
                .repositoryName(repositoryName)
                .folderPath(path ?: "/")
                .commitSpecifier(revision ?: "HEAD")
                .build()
                
            def response = client.getFolder(request)
            
            List<RepositoryEntry> entries = []
            
            // Add files
            response.files()?.each { file ->
                entries.add(new RepositoryEntry(
                    name: file.relativePath().split('/').last(),
                    path: ensureAbsolutePath(file.relativePath()),
                    type: RepositoryProvider.EntryType.FILE,
                    sha: file.blobId(),
                    size: null  // AWS CodeCommit API doesn't provide file size in folder response
                ))
            }
            
            // Add subdirectories - but CodeCommit API has limited support for deep traversal
            response.subFolders()?.each { folder ->
                entries.add(new RepositoryEntry(
                    name: folder.relativePath().split('/').last(),
                    path: ensureAbsolutePath(folder.relativePath()),
                    type: RepositoryProvider.EntryType.DIRECTORY,
                    sha: null, // CodeCommit doesn't provide SHA for directories
                    size: null
                ))
                
                // For recursive listing, we would need additional API calls
                // However, this can be expensive and slow for large repositories
                if (depth != 0 && depth != 1) {
                    try {
                        def subEntries = listDirectory(folder.relativePath(), depth == -1 ? -1 : depth - 1)
                        entries.addAll(subEntries)
                    } catch (Exception e) {
                        // Continue with other directories if one fails
                    }
                }
            }
            
            return entries.sort { it.name }
            
        } catch (Exception e) {
            checkMissingCredsException(e)
            throw new UnsupportedOperationException("Directory listing failed for AWS CodeCommit path: $path - ${e.message}", e)
        }
    }

    protected void checkMissingCredsException(Exception e) {
        final errs = [
                "Failed to connect to service endpoint",
                "Unable to load AWS credentials",
                "The security token included in the request is invalid",
                "The request signature we calculated does not match the signature you provided"]
        if( e !instanceof SdkException )
            return
        if( e instanceof CodeCommitException && e.message?.startsWith("Could not find path") ) {
            // it cannot find the request file
            return
        }
        if( errs.find(it-> e.message?.startsWith(it))) {
            final msg = e.message?.split(/\.|\(|:/)[0].trim()
            throw new MissingCredentialsException("Missing or invalid AWS CodeCommit credentials - $msg", e)
        }
        else {
            throw new AbortOperationException("Unexpected error while connecting repository - $e.message", e)
        }
    }

    /** {@inheritDoc} **/
    // called by AssetManager
    @Override
    void validateRepo() {
        try {
            getRepositoryMetadata()
        }
        catch( IOException e ) {
            throw new AbortOperationException("Cannot access ${getEndpointUrl()} - Make sure a repository exists for it in AWS CodeCommit")
        }
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
        catch (TransportException e) {
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
