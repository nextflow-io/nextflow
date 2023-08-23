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
 *
 */

package nextflow.cli

import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.container.inspect.ContainersInspector

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = "Inspect process settings in a pipeline project")
class CmdInspect extends CmdBase {

    @Override
    String getName() {
        return 'inspect'
    }

    @Parameter(names=['-c','-config'], hidden = true)
    List<String> runConfig

    @Parameter(names=['-format'], description = "Inspect output format. Can be 'json' or 'config'")
    String format = 'json'

    @Parameter(names=['-i','-ignore-errors'], description = 'Ignore errors while inspecting the pipeline')
    boolean ignoreErrors

    @Parameter(names=['-profile'], description = 'Use the given configuration profile(s)')
    String profile

    @Parameter(names=['-w','-await'], description = 'Wait for container images to be available when building images with Wave')
    boolean awaitMode = false

    @Parameter(description = 'Project name or repository url')
    List<String> args

    {
        // enable quiet by default
        CliOptions.quietDefault = true
    }

    @Override
    void run() {
        final target = new CmdRun()
        target.launcher = this.launcher
        target.args = args
        target.profile = this.profile
        target.runConfig = this.runConfig
        target.preview = true
        target.previewAction = this.&applyInspect
        target.ansiLog = false
        // run it
        target.run()
    }

    protected void applyInspect(Session session) {
        // disable wave await mode when running
        if( session.config.wave instanceof Map )
            configAwaitMode(session.config.wave as Map)
        // run the inspector
        new ContainersInspector(session.dag)
                .withFormat(format)
                .withIgnoreErrors(ignoreErrors)
                .printContainers()
    }

    protected void configAwaitMode(Map waveConfig) {
        waveConfig.awaitMode = awaitMode
    }
}
