/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
/**
 * CLI sub-command clone
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Command(name = 'clone', description = "Clone a project into a folder")
class CmdClone extends CmdBase implements HubOptions {

    @Parameters(index = '0', description = 'name of the project to clone')
    String pipeline

    @Parameters(arity = '0..1', description = 'target directory')
    String targetName

    @Option(names = ['-r', '-revision'], description = 'Revision to clone - It can be a git branch, tag or revision number')
    String revision

    @Override
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
        manager.clone(target, revision)
        print "\r"
        println "${manager.project} cloned to: $target"
    }
}
