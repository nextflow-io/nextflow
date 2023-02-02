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
import nextflow.cli.CloneImpl

/**
 * CLI `clone` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Clone a project into a folder')
class CloneCmd extends AbstractCmd implements CloneImpl.Options, HubOptions {

    static public final String NAME = 'clone'

    @Parameter(required = true, description = 'name of the project to clone')
    List<String> args

    @Parameter(names = ['-r'], description = 'Revision to clone - It can be a git branch, tag or revision number')
    String revision

    @Override
    String getName() { NAME }

    @Override
    void run() {
        new CloneImpl(this).run()
    }
}
