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

package nextflow.cli
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.npr.client.RegistryClient
import nextflow.cli.pipeline.PipelineCmdSupport
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import nextflow.scm.HubOptions
import nextflow.util.TestOnly
/**
 * CLI sub-command PULL
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Download or update a project (deprecated -- use `nextflow pipeline install -g`)")
class CmdPull extends CmdBase implements HubAware {

    static final public NAME = 'pull'

    @Parameter(description = 'project name or repository url to pull', arity = 1)
    List<String> args

    @Parameter(names=['-a','-all'], description = 'Update all downloaded projects', arity = 0)
    boolean all

    @Parameter(names=['-r','-revision'], description = 'Revision of the project to pull (either a git branch, tag or commit SHA number)')
    String revision

    @Parameter(names=['-d','-deep'], description = 'Create a shallow clone of the specified depth', hidden = true)
    Integer deep

    @Parameter(names=['-m','-migrate'], description = 'Migrate projects to multi-revision strategy', arity = 0)
    boolean migrate

    @Override
    final String getName() { NAME }

    @TestOnly
    protected File root

    @TestOnly
    protected RegistryClient registryClient

    @Override
    boolean isHidden() { true }

    @Override
    void run() {
        log.warn "The `pull` command is deprecated -- use `nextflow pipeline install -g` instead"

        if( deep )
            log.warn "The `-deep` option is deprecated -- it only applies to projects downloaded with the legacy asset strategy"

        if( !all && !args )
            throw new AbortOperationException('Missing argument')

        if( all && args )
            throw new AbortOperationException('Option `all` requires no arguments')

        if( all && revision )
            throw new AbortOperationException('Option `all` is not compatible with `revision`')

        def list = all ? AssetManager.list() : args.toList()
        if( !list ) {
            log.info "(nothing to do)"
            return
        }

        if( root ) {
            AssetManager.root = root
        }

        // init plugin system
        Plugins.init()

        for( String proj : list ) {
            if( all ) {
                try (def mgr = new AssetManager(proj)) {
                    def branches = mgr.getBranchesAndTags(false).pulled as List<String>
                    branches.each { rev -> pullProjectRevision(proj, rev) }
                }
            } else {
                pullProject(proj, revision)
            }
        }
    }

    private void pullProject(String project, String revision) {
        final address = PipelineCmdSupport.registryAddress(project, hubProvider) { PipelineCmdSupport.registryConfig(this) }
        if( address )
            pullRegistryRelease(address, revision)
        else
            pullProjectRevision(project, revision)
    }

    private void pullRegistryRelease(PipelineCmdSupport.RegistryAddress address, String revision) {
        final release = PipelineCmdSupport.resolveRelease(address, revision, registryClient)
        final result = PipelineCmdSupport.downloadRelease(release, new HubOptions(null, hubUser), deep)
        log.info "${result.message ? ' ' + result.message : ' done'} - revision: ${result.script.revisionInfo}"
    }

    private void pullProjectRevision(String project, String revision) {
        final result = PipelineCmdSupport.downloadRevision(project, revision, toHubOptions(), deep, migrate)
        log.info "${result.message ? ' ' + result.message : ' done'} - revision: ${result.script.revisionInfo}"
    }

}
