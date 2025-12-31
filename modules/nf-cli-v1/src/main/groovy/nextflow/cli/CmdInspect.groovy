/*
 * Copyright 2013-2024, Seqera Labs
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

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileDynamic
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.container.inspect.ContainerInspectMode
import nextflow.container.inspect.ContainersInspector
import nextflow.util.LoggerHelper
/**
 * Implement `inspect` command
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Inspect process settings in a pipeline project")
class CmdInspect extends CmdBase {

    @Override
    String getName() {
        return 'inspect'
    }

    @Parameter(names=['-concretize'], description = "Build the container images resolved by the inspect command")
    boolean concretize

    @Parameter(names=['-c','-config'], hidden = true)
    List<String> runConfig

    @Parameter(names=['-format'], description = "Inspect output format. Can be 'json' or 'config'")
    String format = 'json'

    @Parameter(names=['-i','-ignore-errors'], description = 'Ignore errors while inspecting the pipeline')
    boolean ignoreErrors

    @DynamicParameter(names = '--', hidden = true)
    Map<String,String> params = new LinkedHashMap<>()

    @Parameter(names='-params-file', description = 'Load script parameters from a JSON/YAML file')
    String paramsFile

    @Parameter(names=['-profile'], description = 'Use the given configuration profile(s)')
    String profile

    @Parameter(names=['-r','-revision'], description = 'Revision of the project to inspect (either a git branch, tag or commit SHA number)')
    String revision

    @Parameter(description = 'Project name or repository url')
    List<String> args

    @Override
    void run() {
        ContainerInspectMode.activate(!concretize)
        // configure quiet mode
        LoggerHelper.setQuiet(true)
        // setup the target run command
        final target = new CmdRun()
        target.launcher = this.launcher
        target.args = args
        target.profile = this.profile
        target.revision = this.revision
        target.runConfig = this.runConfig
        target.params = this.params
        target.paramsFile = this.paramsFile
        target.preview = true
        target.previewAction = this.&applyInspect
        target.ansiLog = false
        target.skipHistoryFile = true
        // run it
        target.run()
    }

    protected void applyInspect(Session session) {
        // slow down max rate when concretize is specified
        if( concretize ) {
            configureMaxRate(session.config)
        }
        // run the inspector
        new ContainersInspector(concretize)
                .withFormat(format)
                .withIgnoreErrors(ignoreErrors)
                .printContainers()
    }

    @CompileDynamic
    protected void configureMaxRate(Map config) {
        if( config.wave == null )
            config.wave = new HashMap()
        if( config.wave.httpClient == null )
            config.wave.httpClient = new HashMap()
        config.wave.httpClient.maxRate = '5/30sec'
    }

}
