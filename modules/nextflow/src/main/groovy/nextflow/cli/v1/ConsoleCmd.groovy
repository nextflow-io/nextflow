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
import nextflow.cli.ConsoleImpl

/**
 * CLI `console` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Launch Nextflow interactive console')
class ConsoleCmd extends AbstractCmd implements ConsoleImpl.Options {

    @Parameter(description = 'Nextflow console arguments')
    List<String> args = []

    @Override
    String getScript() {
        args.size() > 0 ? args[0] : null
    }

    @Override
    String getName() { 'console' }

    @Override
    void run() {
        new ConsoleImpl(this).run()
    }
}
