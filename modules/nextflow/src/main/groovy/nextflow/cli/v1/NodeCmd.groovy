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

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.cli.ILauncherOptions
import nextflow.cli.NodeImpl

/**
 * CLI `node` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Launch Nextflow in deamon mode')
class NodeCmd extends AbstractCmd implements NodeImpl.Options {

    @Parameter(names = ['-bg'], arity = 0, description = 'Start the cluster node daemon in background')
    void setBackground(boolean value) {
        launcher.options.background = value
    }

    @DynamicParameter(names = ['-cluster.'], description = 'Define cluster config options')
    Map<String,String> clusterOptions = [:]

    @Parameter(description = 'Daemon name or class')
    List<String> args = []

    @Override
    String getProvider() {
        args.size() ? args[0] : null
    }

    @Override
    ILauncherOptions getLauncherOptions() {
        launcher.options
    }

    @Override
    String getName() { 'node' }

    @Override
    void run() {
        new NodeImpl(this).run()
    }

}
