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

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.container.inspect.ContainersInspector
import nextflow.exception.AbortOperationException
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

    private static final String YES = 'Y'

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

    @Parameter(names=['-r','-revision'], description = 'Revision of the project to inspect (either a git branch, tag or commit SHA number)')
    String revision

    @DynamicParameter(names = '--', hidden = true)
    Map<String,String> params = new LinkedHashMap<>()

    @Parameter(names='-params-file', description = 'Load script parameters from a JSON/YAML file')
    String paramsFile

    @Parameter(names=['-y','-yes'], description = "Automatically reply with 'yes' when prompt for confirmation")
    String assumeYes

    @Parameter(description = 'Project name or repository url')
    List<String> args

    @Override
    void run() {
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
        // run it
        target.run()
    }

    protected void applyInspect(Session session) {
        // disable wave await mode when running
        if( session.config.wave instanceof Map )
            checkWaveConfig(session.config.wave as Map)
        // run the inspector
        new ContainersInspector(session.dag)
                .withFormat(format)
                .withIgnoreErrors(ignoreErrors)
                .printContainers()
    }

    protected void checkWaveConfig(Map wave) {
        if( !wave.enabled || !wave.freeze )
            return 
        if( !assumeYes ) {
            final reply = promptConfirmation()
            if( reply != YES )
                throw new AbortOperationException("Inspect command aborted")
        }
        wave.awaitMode = awaitMode
    }

    protected String promptConfirmation() {
        final console = System.console()
        if( !console ) {
            log.debug "Unable to acquire system console"
            return YES
        }

        print "This command may trigger the build of the container images using Wave. Please confirm the operation [$YES/n]: "
        return console.readLine()
    }
}
