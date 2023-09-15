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
import nextflow.cli.CmdConfig
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import picocli.CommandLine.ParentCommand

/**
 * CLI `config` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'config',
    description = 'Print a project configuration'
)
class ConfigCmd extends AbstractCmd implements CmdConfig.Options {

    @ParentCommand
    private Launcher launcher

    @Parameters(arity = '0..1', description = 'project name')
    String pipeline

    @Option(names = ['-a','--show-profiles'], description = 'Show all configuration profiles')
    boolean showAllProfiles

    @Option(names = ['--profile'], description = 'Choose a configuration profile')
    String profile

    @Option(names = ['--properties'], description = 'Prints config using Java properties notation')
    boolean printProperties

    @Option(names = ['--flat'], description = 'Print config using flat notation')
    boolean printFlatten

    @Option(names = ['--sort'], description = 'Sort config attributes')
    boolean sort

    @Option(names = ['--value'], description = 'Print the value of a config option, or fail if the option is not defined')
    String printValue

    @Override
    CliOptions getLauncherOptions() {
        launcher.options
    }

    @Override
    void run() {
        new CmdConfig(this).run()
    }

}
