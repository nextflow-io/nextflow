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
import nextflow.cli.CmdLog
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters

/**
 * CLI `log` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'log',
    description = 'Print executions log and runtime info'
)
class LogCmd extends AbstractCmd implements CmdLog.Options {

    @Option(names = ['--after'], paramLabel = '<name>|<id>', description = 'Show log entries for runs executed after the specified one')
    String after

    @Option(names = ['--before'], paramLabel = '<name>|<id>', description = 'Show log entries for runs executed before the specified one')
    String before

    @Option(names = ['--but'], paramLabel = '<name>|<id>', description = 'Show log entries of all runs except the specified one')
    String but

    @Option(names = ['-f','--fields'], description = 'Comma separated list of fields to include in the printed log -- Use the `-l` option to show the list of available fields')
    String fields

    @Option(names = ['-F','--filter'], paramLabel = '<expr>', description = "Filter log entries by a custom expression e.g. process =~ /foo.*/ && status == 'COMPLETED'")
    String filterStr

    @Option(names = ['-l','--list-fields'], arity = '0', description = 'Show all available fields')
    boolean listFields

    @Option(names = ['-q','--quiet'], arity = '0', description = 'Show only run names')
    boolean quiet

    @Option(names = ['-s','--separator'], description = 'Character used to separate column values')
    String separator = '\\t'

    @Option(names = ['-t','--template'], paramLabel = '<template>', description = 'Text template used to print each record in the log ')
    String templateStr

    @Parameters(description = 'Session IDs or run names')
    List<String> args = []

    @Override
    void run() {
        new CmdLog(this).run()
    }

}
