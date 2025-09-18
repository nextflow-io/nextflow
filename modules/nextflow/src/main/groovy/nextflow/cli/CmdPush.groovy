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

package nextflow.cli
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.util.TestOnly
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.errors.RepositoryNotFoundException

/**
 * CLI sub-command Push
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Pushes a local implementation to a remote repository")
class CmdPush extends CmdBase implements HubOptions {

    static final public NAME = 'push'

    @Parameter(description = 'Path to push', arity = 1)
    String folderPath

    @Parameter(names=['-repo'], description = 'Defines the repository to push to', required = true)
    String repository

    @Parameter(names=['-r','-revision'], description = 'Revision of the project to run (either a git branch, tag or commit SHA number)')
    String revision = 'main'

    @Parameter(names=['-max-size'], description = 'Maximum file size in MB to push without confirmation (default: 10)')
    int maxSizeMB = 10

    @Parameter(names=['-message', '-m'], description = 'Commit message')
    String message = 'Push from nextflow'

    @Override
    final String getName() { NAME }

    @TestOnly
    protected File root

    @Override
    void run() {
        if( !folderPath )
            throw new AbortOperationException('Missing folder argument')

        def folder = new File(folderPath).getAbsoluteFile()

        if( !folder.exists() )
            throw new AbortOperationException("Folder does not exist: ${folder.absolutePath}")

        if( !folder.isDirectory() )
            throw new AbortOperationException("Path is not a directory: ${folder.absolutePath}")

        log.info "Pushing folder ${folder.absolutePath} to repository ${repository}"

        // init plugin system
        Plugins.init()

        try {
            pushFolder(folder, repository, revision)
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to push folder: ${e.message}", e)
        }
    }

    private void pushFolder(File folder, String repo, String rev) {
        def gitDir = new File(folder, '.git')

        if( gitDir.exists() ) {
            log.debug "Found existing git repository in ${folder.absolutePath}"
            validateExistingRepo(folder, repo)
        } else {
            log.debug "No git repository found, initializing new one"
            initializeRepo(folder, repo, rev)
        }

        checkFileSizes(folder)
        stageAndCommitFiles(folder)
        pushToRemote(folder, rev)

        log.info "Successfully pushed to ${repo} (revision: ${rev})"
    }

    private void validateExistingRepo(File folder, String expectedRepo) {
        try {
            def git = Git.open(folder)
            def config = git.getRepository().getConfig()
            def remoteUrl = config.getString("remote", "origin", "url")

            if( remoteUrl ) {
                def normalizedRemote = normalizeRepoUrl(remoteUrl)
                def normalizedExpected = normalizeRepoUrl(expectedRepo)

                if( normalizedRemote != normalizedExpected ) {
                    throw new AbortOperationException(
                        "Repository mismatch!\n" +
                        "  Local repository: ${remoteUrl}\n" +
                        "  Expected repository: ${expectedRepo}\n" +
                        "Please remove the .git directory or specify the correct repository."
                    )
                }
            }
            git.close()
        }
        catch( RepositoryNotFoundException e ) {
            throw new AbortOperationException("Invalid git repository in ${folder.absolutePath}")
        }
    }

    private String normalizeRepoUrl(String url) {
        return url?.toLowerCase()?.replaceAll(/\.git$/, '')?.replaceAll(/\/$/, '')
    }

    private void initializeRepo(File folder, String repo, String rev) {
        try {
            log.debug "Initializing git repository in ${folder.absolutePath}"
            def git = Git.init().setDirectory(folder).call()

            // Add remote origin
            git.remoteAdd()
                .setName("origin")
                .setUri(new org.eclipse.jgit.transport.URIish(repo))
                .call()

            git.close()
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to initialize git repository: ${e.message}", e)
        }
    }

    private void checkFileSizes(File folder) {
        def maxSizeBytes = maxSizeMB * 1024 * 1024
        List<Map<String,Object>> largeFiles = []

        folder.eachFileRecurse { file ->
            if( file.isFile() && !file.absolutePath.contains('/.git/') ) {
                if( file.length() > maxSizeBytes ) {
                    Map<String,Object> fileEntry = [:]
                    fileEntry.file = file
                    fileEntry.sizeMB = file.length() / (1024 * 1024)
                    largeFiles.add(fileEntry)
                }
            }
        }

        if( largeFiles ) {
            log.warn "Found ${largeFiles.size()} large files:"
            largeFiles.each { entry ->
                def fileInfo = entry.file as File
                def sizeMB = entry.sizeMB as Double
                log.warn "  ${fileInfo.name}: ${String.format('%.1f', sizeMB)} MB"
            }

            print "Do you want to continue and push these large files? [y/N]: "
            def response = System.in.newReader().readLine()?.trim()?.toLowerCase()

            if( response != 'y' && response != 'yes' ) {
                // Add large files to .gitignore
                def fileNames = largeFiles.collect { entry -> (entry.file as File).name }
                addToGitignore(folder, fileNames)
                throw new AbortOperationException("Push cancelled due to large files. Files have been added to .gitignore")
            }
        }
    }

    private void addToGitignore(File folder, List<String> filenames) {
        def gitignoreFile = new File(folder, '.gitignore')
        def content = []

        if( gitignoreFile.exists() ) {
            content = gitignoreFile.readLines()
        }

        filenames.each { filename ->
            if( !content.contains(filename) ) {
                content.add(filename)
            }
        }

        gitignoreFile.text = content.join('\n') + '\n'
        log.info "Added ${filenames.size()} large files to .gitignore"
    }

    private void stageAndCommitFiles(File folder) {
        try {
            def git = Git.open(folder)

            // Add all files
            git.add().addFilepattern(".").call()

            // Check if there are any changes to commit
            def status = git.status().call()
            if( status.clean ) {
                log.info "No changes to commit"
                git.close()
                return
            }

            // Commit changes
            git.commit()
                .setMessage(message)
                .call()

            log.debug "Committed changes with message: ${message}"
            git.close()
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to stage and commit files: ${e.message}", e)
        }
    }

    private void pushToRemote(File folder, String rev) {
        try {
            def git = Git.open(folder)

            // Create and checkout branch if it doesn't exist
            try {
                git.checkout().setName(rev).call()
            }
            catch( Exception ignored ) {
                // Branch doesn't exist, create it
                git.checkout()
                    .setCreateBranch(true)
                    .setName(rev)
                    .call()
            }

            // Push to remote
            def refSpec = "refs/heads/${rev}:refs/heads/${rev}"
            def pushCommand = git.push()
                .setRemote("origin")
                .add(refSpec)

            pushCommand.call()

            log.debug "Push completed successfully"
            git.close()
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to push to remote repository: ${e.message}", e)
        }
    }


}
