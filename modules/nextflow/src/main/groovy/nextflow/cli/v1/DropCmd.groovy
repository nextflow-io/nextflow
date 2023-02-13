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
import nextflow.cli.DropImpl

/**
 * CLI `drop` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Delete the local copy of a project')
class DropCmd extends AbstractCmd implements DropImpl.Options {

    @Parameter(required = true, description = 'name of the project to drop')
    List<String> args

    @Parameter(names = ['-f','-force'], description = 'Delete the repository without taking care of local changes')
    boolean force

    @Override
    String getPipeline() { args[0] }

    @Override
    String getName() { 'drop' }

    @Override
    void run() {
        new DropImpl(this).run()
    }
}
