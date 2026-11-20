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

package nextflow.cli.project

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.cli.CmdBase
import nextflow.plugin.Plugins
import nextflow.scm.AssetManager

/**
 * Project remove subcommand -- delete the local copy of a project
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Delete the local copy of a project")
class CmdProjectRemove extends CmdBase {

    static final public String NAME = 'remove'

    @Parameter(required=true, description = 'name of the project to remove')
    List<String> args

    @Parameter(names=['-r','-revision'], description = 'Revision of the project to remove (either a git branch, tag or commit SHA number)')
    String revision

    @Parameter(names=['-f','-force'], description = 'Delete the repository without taking care of local changes', arity = 0)
    boolean force

    @Override
    String getName() { NAME }

    @Override
    void run() {
        Plugins.init()
        try (final manager = new AssetManager(args[0])) {
            manager.drop(revision, force)
        }
    }
}
