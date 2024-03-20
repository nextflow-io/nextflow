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
import nextflow.cli.CmdView
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * CLI `view` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'view',
    description = 'View project script file(s)'
)
class ViewCmd extends AbstractCmd implements CmdView.Options {

    @Parameters(description = 'project name')
    String pipeline

    @Option(names = ['-l','--all'], arity = '0', description = 'List repository content')
    boolean all

    @Option(names = ['-q','--quiet'], arity = '0', description = 'Hide header line')
    boolean quiet

    @Override
    void run() {
        new CmdView(this).run()
    }
}
