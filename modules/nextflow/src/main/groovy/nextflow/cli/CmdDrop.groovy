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
 * CLI sub-command DROP
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Delete the local copy of a project")
class CmdDrop extends CmdBase {

    static final public NAME = 'drop'

    @Parameter(required=true, description = 'name of the project to drop')
    List<String> args

    @Parameter(names='-f', description = 'Delete the repository without taking care of local changes')
    boolean force

    @Override
    final String getName() { NAME }

    @Override
    void run() {

        def manager = new AssetManager(args[0])
        if( !manager.localPath.exists() ) {
            throw new AbortOperationException("No match found for: ${args[0]}")
        }

        if( this.force || manager.isClean() ) {
            manager.close()
            if( !manager.localPath.deleteDir() )
                throw new AbortOperationException("Unable to delete project `${manager.project}` -- Check access permissions for path: ${manager.localPath}")
            return
        }

        throw new AbortOperationException("Local project repository contains uncommitted changes -- wont drop it")
    }
}
