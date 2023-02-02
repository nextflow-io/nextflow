/*
 * Copyright 2023, Seqera Labs
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
import nextflow.cli.CloneImpl
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
class CloneCmd extends AbstractCmd implements CloneImpl.Options, HubOptions {

    @Parameters(arity = '1..', description = 'name of the project to clone')
    List<String> args

    @Option(names = ['-r','--revision'], description = 'Revision to clone - It can be a git branch, tag or revision number')
    String revision

    @Override
    Integer call() {
        new CloneImpl(this).run()
        return 0
    }
}
