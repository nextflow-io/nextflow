/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.cli.v1

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.cli.CleanImpl
import nextflow.cli.ILauncherOptions

/**
 * CLI `clean` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Lorenz Gerber <lorenzottogerber@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Clean up project cache and work directories')
class CleanCmd extends AbstractCmd implements CleanImpl.Options {

    static public final String NAME = 'clean'

    @Parameter(names = ['-after'], description = 'Clean up runs executed after the specified one')
    String after

    @Parameter(names = ['-before'], description = 'Clean up runs executed before the specified one')
    String before

    @Parameter(names = ['-but'], description = 'Clean up all runs except the specified one')
    String but

    @Parameter(names = ['-n', '-dry-run'], arity = 0, description = 'Print names of file to be removed without deleting them')
    boolean dryRun

    @Parameter(names = ['-f', '-force'], arity = 0, description = 'Force clean command')
    boolean force

    @Parameter(names = ['-k', '-keep-logs'], description = 'Removes only temporary files but retains execution log entries and metadata')
    boolean keepLogs

    @Parameter(names = ['-q', '-quiet'], arity = 0, description = 'Do not print names of files removed')
    boolean quiet

    @Parameter
    List<String> args

    @Override
    ILauncherOptions getLauncherOptions() {
        launcher.options
    }

    @Override
    String getName() { NAME }

    @Override
    void run() {
        new CleanImpl(this).run()
    }

}
