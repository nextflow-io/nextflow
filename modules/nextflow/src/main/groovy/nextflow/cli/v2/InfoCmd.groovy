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
import nextflow.cli.CmdInfo
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * CLI `info` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'info',
    description = 'Print project and system runtime information'
)
class InfoCmd extends AbstractCmd implements CmdInfo.Options {

    @Parameters(arity = '0..1', description = 'project name')
    String pipeline

    @Option(names = ['-d'], arity = '0', description = 'Show detailed information')
    boolean detailed

    @Option(names = ['-dd'], arity = '0', hidden = true)
    boolean moreDetailed

    @Option(names = ['-o','--output-format'], description = 'Output format, either: text (default), json, yaml')
    String format

    @Option(names = ['-u','--check-updates'], description = 'Check for remote updates')
    boolean checkForUpdates

    @Override
    void run() {
        new CmdInfo(this).run()
    }

}
