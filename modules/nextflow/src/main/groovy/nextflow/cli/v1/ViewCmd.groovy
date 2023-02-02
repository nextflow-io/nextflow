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
import nextflow.cli.ViewImpl

/**
 * CLI `view` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'View project script file(s)')
class ViewCmd extends AbstractCmd implements ViewImpl.Options {

    static public final String NAME = 'view'

    @Override
    String getName() { NAME }

    @Parameter(description = 'project name', required = true)
    List<String> args = []

    @Parameter(names = ['-l'], arity = 0, description = 'List repository content')
    boolean all

    @Parameter(names = ['-q'], arity = 0, description = 'Hide header line')
    boolean quiet

    @Override
    void run() {
        new ViewImpl(this).run()
    }
}
