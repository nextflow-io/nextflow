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
import nextflow.cli.CmdPull
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * CLI `pull` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'pull',
    description = 'Download or update a project'
)
class PullCmd extends AbstractCmd implements CmdPull.Options, HubOptionsV2 {

    @Parameters(description = 'project name or repository url to pull')
    String pipeline

    @Option(names = ['--all'], arity = '0', description = 'Update all downloaded projects')
    boolean all

    @Option(names = ['-d','--depth'], description = 'Create a shallow clone of the specified depth')
    Integer depth

    @Option(names = ['-r','--revision'], description = 'Revision of the project to run (either a git branch, tag or commit SHA number)')
    String revision

    @Override
    void run() {
        new CmdPull(this).run()
    }

}
