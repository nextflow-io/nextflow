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
import nextflow.cli.CmdNode
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import picocli.CommandLine.ParentCommand

/**
 * CLI `node` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'node',
    description = 'Launch Nextflow in daemon mode'
)
class NodeCmd extends AbstractCmd implements CmdNode.Options {

    @ParentCommand
    private Launcher launcher

    @Option(names = ['--cluster.'], paramLabel = '<name>=<value>', description = 'Define cluster config options')
    Map<String,String> clusterOptions = [:]

    @Option(names = ['--bg'], arity = '0', description = 'Start the cluster node daemon in background')
    void setBackground(boolean value) {
        launcher.options.background = value
    }

    @Parameters(description = 'Daemon name or class')
    String provider

    @Override
    CliOptions getLauncherOptions() {
        launcher.options
    }

    @Override
    void run() {
        new CmdNode(this).run()
    }

}
