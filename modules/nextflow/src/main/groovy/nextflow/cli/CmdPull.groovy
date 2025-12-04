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
import nextflow.scm.MultiRevisionAssetManager
import nextflow.util.TestOnly
/**
 * CLI sub-command PULL
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Download or update a project")
class CmdPull extends CmdBase implements HubOptions {

    static final public NAME = 'pull'

    @Parameter(description = 'project name or repository url to pull', arity = 1)
    List<String> args

    @Parameter(names=['-a','-all'], description = 'Update all downloaded projects', arity = 0)
    boolean all

    @Parameter(names=['-r','-revision'], description = 'Revision of the project to pull (either a git branch, tag or commit SHA number)')
    String revision

    @Override
    final String getName() { NAME }

    @TestOnly
    protected File root

    @Override
    void run() {

        if( !all && !args )
            throw new AbortOperationException('Missing argument')

        if( all && args )
            throw new AbortOperationException('Option `all` requires no arguments')

        if( all && revision )
            throw new AbortOperationException('Option `all` is not compatible with `revision`')

        if( root ) {
            AssetManager.root = root
        }

        // init plugin system
        Plugins.init()

        List<MultiRevisionAssetManager> list = []
        if ( all ) {
            def all = AssetManager.list()
            all.each{ proj ->
                def revManager = new MultiRevisionAssetManager(proj)
                revManager.listRevisions().each{ rev ->
                    list << new MultiRevisionAssetManager(proj, this, rev)
                }
            }
        } else {
            args.toList().each {
                list << new MultiRevisionAssetManager(it, this, revision)
            }
        }

        if( !list ) {
            log.info "(nothing to do)"
            return
        }

        list.each { manager ->
            log.info "Checking ${manager.getProjectWithRevision()} ..."

            def result = manager.createSharedClone()
            manager.updateModules()

            def scriptFile = manager.getScriptFile()
            String message = !result ? " done" : " $result"
            message += " - revision: ${scriptFile.revisionInfo}"
            log.info message
        }

    }

}
