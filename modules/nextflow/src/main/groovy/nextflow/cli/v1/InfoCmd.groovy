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
import nextflow.cli.InfoImpl

/**
 * CLI `info` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Print project and system runtime information')
class InfoCmd extends AbstractCmd implements InfoImpl.Options {

    static public final String NAME = 'info'

    @Parameter(description = 'project name')
    List<String> args

    @Parameter(names = ['-d'], arity = 0, description = 'Show detailed information')
    boolean detailed

    @Parameter(names = ['-dd'], arity = 0, hidden = true)
    boolean moreDetailed

    @Parameter(names = ['-o'], description = 'Output format, either: text (default), json, yaml')
    String format

    @Parameter(names = ['-u','-check-updates'], description = 'Check for remote updates')
    boolean checkForUpdates

    @Override
    String getName() { NAME }

    @Override
    void run() {
        new InfoImpl(this).run()
    }

}
