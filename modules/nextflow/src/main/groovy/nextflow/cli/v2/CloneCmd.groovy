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

package nextflow.cli.v2

import groovy.transform.CompileStatic
import nextflow.cli.CmdClone
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * CLI `clone` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'clone',
    description = 'Clone a project into a folder'
)
class CloneCmd extends AbstractCmd implements CmdClone.Options, HubOptionsV2 {

    @Parameters(index = '0', description = 'name of the project to clone')
    String pipeline

    @Parameters(arity = '0..1', paramLabel = '<target>', description = 'target directory')
    String targetName

    @Option(names = ['-d','--depth'], description = 'Create a shallow clone of the specified depth')
    Integer depth

    @Option(names = ['-r','--revision'], description = 'Revision to clone - It can be a git branch, tag or revision number')
    String revision

    @Override
    void run() {
        new CmdClone(this).run()
    }
}
