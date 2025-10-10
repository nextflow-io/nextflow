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
import nextflow.scm.AssetManager
import nextflow.util.TestOnly
import org.eclipse.jgit.api.Git
import org.eclipse.jgit.transport.RemoteConfig
import org.eclipse.jgit.transport.URIish


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

    @Parameter(description = 'Repository URL to push to (optional if already configured as git remote)')
    List<String> args

    @Parameter(names=['-d', '-directory'], description = 'Local directory to push (default: current directory)')
    String directory

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
        if( args && args.size() > 1){
            throw new AbortOperationException('Incorrect number of arguments')
        }

        // Get repository from args (optional)
        def repository = args && args.size() == 1 ? args[0] : null

        // Folder defaults to current working directory if not specified
        def folder = directory
            ? new File(directory).getAbsoluteFile()
            : new File(System.getProperty('user.dir')).getAbsoluteFile()

        if( !folder.exists() )
            throw new AbortOperationException("Folder does not exist: ${folder.absolutePath}")

        if( !folder.isDirectory() )
            throw new AbortOperationException("Path is not a directory: ${folder.absolutePath}")

        // init plugin system
        Plugins.init()

        try {
            def resolvedRepo = repository
            if( !resolvedRepo ) {
                resolvedRepo = resolveRepository(folder)
            }

            log.info "Pushing folder ${folder.absolutePath} to repository ${resolvedRepo}"
            pushFolder(folder, resolvedRepo, revision)
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to push folder: ${e.message}", e)
        }
    }

    private void pushFolder(File folder, String repo, String rev) {
        def gitDir = new File(folder, '.git')
        def remoteName = "origin"
        def isNewRepo = false

        if( gitDir.exists() ) {
            log.debug "Found existing git repository in ${folder.absolutePath}"
            remoteName = validateExistingRepo(folder, repo)
            checkCurrentBranch(folder, rev)
        } else {
            log.debug "No git repository found, initializing new one"
            initializeRepo(folder, repo, rev)
            isNewRepo = true
        }

        checkFileSizes(folder)
        manageNextflowGitignore(folder)
        stageAndCommitFiles(folder)
        def manager = new AssetManager(folder, repo, this)
        manager.upload(rev, remoteName, isNewRepo)
        log.info "Successfully pushed to ${repo} (revision: ${rev})"
    }

    private String validateExistingRepo(File folder, String expectedRepo) {
        def git = Git.open(folder)

        try {
            def remotes = git.remoteList().call()

            // Find all remotes and check if any matches the expected repo
            def matchingRemote = null

            for( RemoteConfig remote : remotes ) {
                if( remote.URIs ) {
                    def remoteUrl = remote.URIs[0].toString()
                    def normalizedRemote = normalizeRepoUrl(remoteUrl)
                    def normalizedExpected = normalizeRepoUrl(expectedRepo)

                    if( normalizedRemote == normalizedExpected ) {
                        matchingRemote = remote.name
                        break
                    }
                }
            }

            if( !matchingRemote ) {
                def remotesList = remotes.collect { remote ->
                    def url = remote.URIs ? remote.URIs[0].toString() : 'no URL'
                    "  ${remote.name}: ${url}"
                }.join('\n')

                throw new AbortOperationException(
                    "Repository URL not found in remotes!\n" +
                    "  Expected repository: ${expectedRepo}\n" +
                    "  Available remotes:\n${remotesList}\n" +
                    "Please add the repository as a remote or specify the correct repository."
                )
            }

            return matchingRemote
        }
        finally {
            git.close()
        }
    }

    private String normalizeRepoUrl(String url) {
        return url?.toLowerCase()?.replaceAll(/\.git$/, '')?.replaceAll(/\/$/, '')
    }

    private void checkCurrentBranch(File folder, String requestedBranch) {
        def git = Git.open(folder)

        try {
            def head = git.getRepository().findRef("HEAD")
            if( !head ) {
                log.debug "No HEAD found, assuming new repository"
                git.close()
                return
            }

            def currentBranch = null
            if( head.isSymbolic() ) {
                currentBranch = git.getRepository().getBranch()
            } else {
                log.debug "HEAD is not symbolic (detached state)"
                git.close()
                throw new AbortOperationException("Repository is in detached HEAD state. Please checkout to a branch before pushing.")
            }

            if( currentBranch && currentBranch != requestedBranch ) {
                git.close()
                throw new AbortOperationException(
                    "Current branch '${currentBranch}' does not match requested branch '${requestedBranch}'.\n" +
                    "Please checkout to branch '${requestedBranch}' before pushing or specify the correct branch with -r option."
                )
            }

            log.debug "Current branch '${currentBranch}' matches requested branch '${requestedBranch}'"
        }
        finally {
            git.close()
        }
    }

    private void initializeRepo(File folder, String repo, String rev) {
        log.debug "Initializing git repository in ${folder.absolutePath}"
        def git = Git.init().setDirectory(folder).call()

        // Add remote origin
        git.remoteAdd()
            .setName("origin")
            .setUri(new URIish(repo))
            .call()

        git.close()
    }

    private void checkFileSizes(File folder) {
        def maxSizeBytes = maxSizeMB * 1024 * 1024
        def git = Git.open(folder)

        try {
            // Get Git status to find files that would be committed
            def status = git.status().call()
            def filesToBeCommitted = []

            // Add untracked files
            filesToBeCommitted.addAll(status.untracked)
            // Add modified files
            filesToBeCommitted.addAll(status.modified)
            // Add added files
            filesToBeCommitted.addAll(status.added)

            def largeFiles = []

            filesToBeCommitted.each { relativePath ->
                def file = new File(folder, relativePath as String)
                if( file.exists() && file.isFile() && file.length() > maxSizeBytes ) {
                    def fileEntry = [
                        file: file,
                        relativePath: relativePath,
                        sizeMB: file.length() / (1024 * 1024)
                    ]
                    largeFiles.add(fileEntry)
                }
            }

            if( largeFiles ) {
                log.warn "Found ${largeFiles.size()} large files that would be committed:"
                largeFiles.each { entry ->
                    def sizeMB = entry['sizeMB'] as Double
                    log.warn "  ${entry['relativePath']}: ${String.format('%.1f', sizeMB)} MB"
                }

                print "Do you want to push these large files? [y/N]: "
                def response = System.in.newReader().readLine()?.trim()?.toLowerCase()

                if( response != 'y' && response != 'yes' ) {
                    // Add large files to .gitignore
                    def relativePaths = largeFiles.collect { entry -> entry['relativePath'] as String }
                    addToGitignore(folder, relativePaths)
                    println "Files have been added to .gitignore"
                }
            }
        }
        finally {
            git.close()
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

    private void manageNextflowGitignore(File folder) {
        def gitignoreFile = new File(folder, '.gitignore')
        List<String> content = []

        if( gitignoreFile.exists() ) {
            content = gitignoreFile.readLines()
        }

        // Default Nextflow entries to add
        def nextflowEntries = [
            '.nextflow',
            '.nextflow.log*'
        ]

        def added = []
        nextflowEntries.each { entry ->
            if( !content.contains(entry) ) {
                content.add(entry)
                added.add(entry)
            }
        }

        // Check for work directory
        def workDirs = findWorkDirectories(folder)
        if( workDirs ) {
            def workEntriesToAdd = promptForWorkDirectories(workDirs, content)
            workEntriesToAdd.each { workDir ->
                if( !content.contains(workDir) ) {
                    content.add(workDir)
                    added.add(workDir)
                }
            }
        }

        if( added ) {
            gitignoreFile.text = content.join('\n') + '\n'
            log.info "Added ${added.size()} Nextflow entries to .gitignore: ${added.join(', ')}"
        } else {
            log.debug "All Nextflow entries already present in .gitignore"
        }
    }

    private List<String> findWorkDirectories(File folder) {
        List<String> workDirs = []

        // Check for the default Nextflow work directory
        def workDir = new File(folder, 'work')
        if( workDir.exists() && workDir.isDirectory() ) {
            workDirs.add('work')
        }

        return workDirs
    }

    private List<String> promptForWorkDirectories(List<String> workDirs, List<String> currentGitignore) {
        List<String> toAdd = []

        workDirs.each { workDir ->
            // Check if already in .gitignore
            if( currentGitignore.contains(workDir) ) {
                log.debug "Work directory '${workDir}' already in .gitignore"
                return // Skip this directory
            }

            println "Found Nextflow work directory: ${workDir}"
            print "Do you want to add '${workDir}' to .gitignore? [Y/n]: "
            def response = System.in.newReader().readLine()?.trim()?.toLowerCase()

            // Default to 'yes' if empty response or 'y'/'yes'
            if( !response || response == 'y' || response == 'yes' ) {
                toAdd.add(workDir)
                log.info "Will add '${workDir}' to .gitignore"
            } else {
                log.info "Skipping '${workDir}'"
            }
        }

        return toAdd
    }

    private void stageAndCommitFiles(File folder) {
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

    private String resolveRepository(File folder) {
        def gitDir = new File(folder, '.git')

        if( !gitDir.exists() ) {
            throw new AbortOperationException("No git repository found and no repository URL provided. Please specify a repository with -repo parameter.")
        }

        def git = Git.open(folder)

        try {
            def remotes = git.remoteList().call()

            if( remotes.empty ) {
                throw new AbortOperationException("No remotes configured in git repository. Please add a remote or specify a repository with -repo parameter.")
            }

            if( remotes.size() == 1 ) {
                def remote = remotes[0]
                def remoteUrl = remote.URIs[0].toString()
                log.info "Using remote '${remote.name}': ${remoteUrl}"
                return remoteUrl
            }

            // Multiple remotes - ask user to choose
            return selectRemoteFromUser(remotes)
        }
        finally {
            git.close()
        }
    }

    private static String selectRemoteFromUser(List<RemoteConfig> remotes) {
        println "Multiple remotes found. Please select which remote to push to:"

        def remoteOptions = [:]
        remotes.eachWithIndex { remote, index ->
            def remoteUrl = remote.URIs[0].toString()
            def remoteInfo = [name: remote.name, url: remoteUrl]
            remoteOptions[index + 1] = remoteInfo
            println "  ${index + 1}. ${remote.name}: ${remoteUrl}"
        }

        println "  ${remotes.size() + 1}. Cancel"

        while( true ) {
            print "Enter your choice [1-${remotes.size() + 1}]: "
            def input = System.in.newReader().readLine()?.trim()

            try {
                def choice = Integer.parseInt(input)

                if( choice == remotes.size() + 1 ) {
                    throw new AbortOperationException("Push operation cancelled by user.")
                }

                if( choice >= 1 && choice <= remotes.size() ) {
                    def selected = remoteOptions[choice]
                    log.info "Selected remote '${selected['name']}': ${selected['url']}"
                    return selected['url']
                }

                println "Invalid choice. Please enter a number between 1 and ${remotes.size() + 1}."
            }
            catch( NumberFormatException ignored ) {
                println "Invalid input. Please enter a number."
            }
        }
    }


}
