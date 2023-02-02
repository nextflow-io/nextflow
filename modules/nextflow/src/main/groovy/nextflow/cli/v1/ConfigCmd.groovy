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
import nextflow.cli.ConfigImpl
import nextflow.cli.ILauncherOptions

/**
 * CLI `config` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Print a project configuration')
class ConfigCmd extends AbstractCmd implements ConfigImpl.Options {

    static public final String NAME = 'config'

    @Parameter(description = 'project name')
    List<String> args = []

    @Parameter(names = ['-a','-show-profiles'], description = 'Show all configuration profiles')
    boolean showAllProfiles

    @Parameter(names = ['-profile'], description = 'Choose a configuration profile')
    String profile

    @Parameter(names = ['-properties'], description = 'Prints config using Java properties notation')
    boolean printProperties

    @Parameter(names = ['-flat'], description = 'Print config using flat notation')
    boolean printFlatten

    @Parameter(names = ['-sort'], description = 'Sort config attributes')
    boolean sort

    @Override
    ILauncherOptions getLauncherOptions() {
        launcher.options
    }

    @Override
    String getName() { NAME }

    @Override
    void run() {
        new ConfigImpl(this).run()
    }

}
