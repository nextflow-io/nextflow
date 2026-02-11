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
import nextflow.scm.PushManager
import nextflow.util.TestOnly


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

    @Parameter(names=['-c', '-commit'], description = 'Add and commit changes in the current directory (default false)')
    boolean commit = false

    @Parameter(names=['-r','-revision'], description = 'Revision of the project to run (either a git branch, tag or commit SHA number)')
    String revision

    @Parameter(names=['-max-size'], description = 'Maximum file size in MB to push without confirmation (default: 10)')
    int maxSizeMB = 10

    @Parameter(names=['-message', '-m'], description = 'Commit message')
    String message = 'Push from nextflow'

    @Parameter(names=['-y', '-yes'], description = "Automatically reply 'yes' to commit prompts")
    boolean autoAccept = false

    @Override
    final String getName() { NAME }

    @TestOnly
    protected File rootFolder

    @Override
    void run() {
        if( args && args.size() > 1){
            throw new AbortOperationException('Incorrect number of arguments')
        }

        // Get repository from args (optional)
        def repository = args && args.size() == 1 ? args[0] : null
        def folder = rootFolder ?: new File(System.getProperty('user.dir')).getAbsoluteFile()

        if( !folder.exists() )
            throw new AbortOperationException("Folder does not exist: ${folder.absolutePath}")

        if( !folder.isDirectory() )
            throw new AbortOperationException("Path is not a directory: ${folder.absolutePath}")

        // init plugin system
        Plugins.init()

        try {
            final manager = new PushManager(folder, commit, maxSizeMB, autoAccept, true)
            def resolvedRepo = repository
            if( !resolvedRepo ) {
                resolvedRepo = manager.resolveRepository()
            }

            log.info "Pushing folder ${folder.absolutePath} to repository ${resolvedRepo} (commit: $commit)"
            manager.push(resolvedRepo, revision)
        }
        catch( Exception e ) {
            throw new AbortOperationException("Failed to push folder: ${e.message}", e)
        }
    }

}
