/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
import nextflow.scm.AssetManager
/**
 * CLI sub-command clone
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Clone a project into a folder")
class CmdClone extends CmdBase implements HubOptions {

    static final public NAME = 'clone'

    @Parameter(required=true, description = 'name of the project to clone')
    List<String> args

    @Parameter(names='-r', description = 'Revision to clone - It can be a git branch, tag or revision number')
    String revision

    @Override
    final String getName() { NAME }

    @Override
    void run() {
        // the pipeline name
        String pipeline = args[0]
        final manager = new AssetManager(pipeline, this)

        // the target directory is the second parameter
        // otherwise default the current pipeline name
        def target = new File(args.size()> 1 ? args[1] : manager.getBaseName())
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
