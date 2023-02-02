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
import nextflow.cli.LogImpl

/**
 * CLI `log` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Print executions log and runtime info')
class LogCmd extends AbstractCmd implements LogImpl.Options {

    static public final String NAME = 'log'

    @Parameter(names = ['-after'], description = 'Show log entries for runs executed after the specified one')
    String after

    @Parameter(names = ['-before'], description = 'Show log entries for runs executed before the specified one')
    String before

    @Parameter(names = ['-but'], description = 'Show log entries of all runs except the specified one')
    String but

    @Parameter(names = ['-f','-fields'], description = 'Comma separated list of fields to include in the printed log -- Use the `-l` option to show the list of available fields')
    String fields

    @Parameter(names = ['-F','-filter'], description = "Filter log entries by a custom expression e.g. process =~ /foo.*/ && status == 'COMPLETED'")
    String filterStr

    @Parameter(names = ['-l','-list-fields'], arity = 0, description = 'Show all available fields')
    boolean listFields

    @Parameter(names = ['-q','-quiet'], arity = 0, description = 'Show only run names')
    boolean quiet

    @Parameter(names = ['-s'], description = 'Character used to separate column values')
    String separator = '\\t'

    @Parameter(names = ['-t','-template'], description = 'Text template used to each record in the log ')
    String templateStr

    @Parameter(description = 'Run name or session id')
    List<String> args

    @Override
    String getName() { NAME }

    @Override
    void run() {
        new LogImpl(this).run()
    }

}
