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
import nextflow.cli.CliOptions
import nextflow.cli.CmdClean
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import picocli.CommandLine.ParentCommand

/**
 * CLI `clean` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'clean',
    description = 'Clean up project cache and work directories'
)
class CleanCmd extends AbstractCmd implements CmdClean.Options {

    @ParentCommand
    private Launcher launcher

    @Option(names = ['--after'], paramLabel = '<name>|<id>', description = 'Clean up runs executed after the specified one')
    String after

    @Option(names = ['--before'], paramLabel = '<name>|<id>', description = 'Clean up runs executed before the specified one')
    String before

    @Option(names = ['--but'], paramLabel = '<name>|<id>', description = 'Clean up all runs except the specified one')
    String but

    @Option(names = ['-n', '--dry-run'], arity = '0', description = 'Print names of file to be removed without deleting them')
    boolean dryRun

    @Option(names = ['-f', '--force'], arity = '0', description = 'Force clean command')
    boolean force

    @Option(names = ['-k', '--keep-logs'], description = 'Removes only temporary files but retains execution log entries and metadata')
    boolean keepLogs

    @Option(names = ['-q', '--quiet'], arity = '0', description = 'Do not print names of files removed')
    boolean quiet

    @Parameters(description = 'Session IDs or run names')
    List<String> args = []

    @Override
    CliOptions getLauncherOptions() {
        launcher.options
    }

    @Override
    void run() {
        new CmdClean(this).run()
    }

}
