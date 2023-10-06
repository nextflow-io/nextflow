/*
 * Copyright 2013-2023, Seqera Labs
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
/**
 * CLI sub-command clone
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class CmdClone {

    static final public NAME = 'clone'

    interface Options extends HubOptions {
        String getPipeline()
        String getTargetName()
        Integer getDepth()
        String getRevision()
    }

    @Parameters(commandDescription = "Clone a project into a folder")
    static class V1 extends CmdBase implements Options, HubOptions.V1 {

        @Parameter(required=true, description = 'name of the project to clone')
        List<String> args = []

        @Parameter(names='-r', description = 'Revision to clone - It can be a git branch, tag or revision number')
        String revision

        @Parameter(names=['-d','-depth','-deep'], description = 'Create a shallow clone of the specified depth')
        Integer depth

        @Override
        String getPipeline() { args[0] }

        @Override
        String getTargetName() {
            args.size() > 1 ? args[1] : null
        }

        @Override
        final String getName() { NAME }

        @Override
        void run() {
            new CmdClone(this).run()
        }
    }

    @Delegate
    private Options options

    CmdClone(Options options) {
        this.options = options
    }

    void run() {
        // init plugin system
        Plugins.init()
        final manager = new AssetManager(pipeline, this)

        // the target directory is the second parameter
        // otherwise default the current pipeline name
        def target = new File(targetName ?: manager.getBaseName())
        if( target.exists() ) {
            if( target.isFile() )
                throw new AbortOperationException("A file with the same name already exists: $target")
            if( !target.empty() )
                throw new AbortOperationException("Clone target directory must be empty: $target")
        }
        else if( !target.mkdirs() ) {
            throw new AbortOperationException("Cannot create clone target directory: $target")
        }

        manager.checkValidRemoteRepo()
        print "Cloning ${manager.project}${revision ? ':'+revision:''} ..."
        manager.clone(target, revision, depth)
        print "\r"
        println "${manager.project} cloned to: $target"
    }
}
