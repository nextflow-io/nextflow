/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.cli

import com.beust.jcommander.DynamicParameter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import org.pf4j.ExtensionPoint

/**
 * Implements the 'nextflow launch' command
 *
 * @author Phil Ewels <phil.ewels@seqera.io>
 */
@Slf4j
@CompileStatic
@Parameters(commandDescription = "Launch a workflow in Seqera Platform")
class CmdLaunch extends CmdBase implements UsageAware {

    interface LaunchCommand extends ExtensionPoint {
        void launch(LaunchOptions options)
    }

    static final public String NAME = 'launch'

    @Override
    String getName() { NAME }

    @Parameter(description = 'Pipeline repository URL')
    List<String> args

    @Parameter(names = ['-workspace'], description = 'Workspace name')
    String workspace

    @Parameter(names = ['-compute-env'], description = 'Compute environment name')
    String computeEnv

    @Parameter(names = ['-name'], description = 'Assign a mnemonic name to the pipeline run')
    String runName

    @Parameter(names = ['-w', '-work-dir'], description = 'Directory where intermediate result files are stored')
    String workDir

    @Parameter(names = ['-r', '-revision'], description = 'Revision of the project to run (either a git branch, tag or commit SHA number)')
    String revision

    @Parameter(names = ['-profile'], description = 'Choose a configuration profile')
    String profile

    @Parameter(names = ['-c', '-config'], description = 'Add the specified file to configuration set')
    List<String> configFiles

    @Parameter(names = ['-params-file'], description = 'Load script parameters from a JSON/YAML file')
    String paramsFile

    @Parameter(names = ['-entry'], description = 'Entry workflow name to be executed')
    String entryName

    @Parameter(names = ['-resume'], description = 'Execute the script using the cached results')
    String resume

    @Parameter(names = ['-latest'], description = 'Pull latest changes before run')
    boolean latest

    @Parameter(names = ['-stub-run', '-stub'], description = 'Execute the workflow replacing process scripts with command stubs')
    boolean stubRun

    @Parameter(names = ['-main-script'], description = 'The script file to be executed when launching a project directory or repository')
    String mainScript

    /**
     * Defines the parameters to be passed to the pipeline script
     */
    @DynamicParameter(names = '--', description = 'Set a parameter used by the pipeline', hidden = true)
    Map<String, String> params = new LinkedHashMap<>()

    private LaunchCommand operation

    @Override
    void run() {
        // Validate required parameters
        if (!args || args.isEmpty()) {
            log.debug "No pipeline repository URL provided"
            throw new AbortOperationException("Pipeline repository URL is required")
        }

        // Load the Launch command implementation
        this.operation = loadOperation()
        if (!operation)
            throw new IllegalStateException("Unable to load launch extension.")

        // Create options object
        final options = new LaunchOptions(
            pipeline: args[0],
            workspace: workspace,
            computeEnv: computeEnv,
            runName: runName,
            workDir: workDir,
            revision: revision,
            profile: profile,
            configFiles: configFiles,
            paramsFile: paramsFile,
            entryName: entryName,
            resume: resume,
            latest: latest,
            stubRun: stubRun,
            mainScript: mainScript,
            params: params,
            launcher: launcher
        )

        // Execute launch
        operation.launch(options)
    }

    @Override
    void usage() {
        usage(null)
    }

    @Override
    void usage(List<String> args) {
        def result = []
        result << this.getClass().getAnnotation(Parameters).commandDescription()
        result << 'Usage: nextflow launch <pipeline> [options]'
        result << ''
        result << 'Options:'
        result << '  -workspace <name>         Workspace name'
        result << '  -compute-env <name>       Compute environment name (default: primary)'
        result << '  -name <name>              Assign a mnemonic name to the pipeline run'
        result << '  -w, -work-dir <path>      Directory where intermediate result files are stored'
        result << '  -r, -revision <revision>  Revision of the project to run (git branch, tag or commit SHA)'
        result << '  -profile <profile>        Choose a configuration profile'
        result << '  -c, -config <file>        Add the specified file to configuration set'
        result << '  -params-file <file>       Load script parameters from a JSON/YAML file'
        result << '  -entry <name>             Entry workflow name to be executed'
        result << '  -resume [session]         Execute the script using the cached results'
        result << '  -latest                   Pull latest changes before run'
        result << '  -stub-run, -stub          Execute the workflow replacing process scripts with command stubs'
        result << '  -main-script <file>       The script file to be executed when launching a project'
        result << '  --<param>=<value>         Set a parameter used by the pipeline'
        result << ''
        println result.join('\n').toString()
    }

    protected LaunchCommand loadOperation() {
        // Setup the plugins system and load the launch provider
        Plugins.init()
        // Load the config
        Plugins.start('nf-tower')
        // Get Launch command operations implementation from plugins
        return Plugins.getExtension(LaunchCommand)
    }

    /**
     * Data class to hold launch options
     */
    @CompileStatic
    static class LaunchOptions {
        String pipeline
        String workspace
        String computeEnv
        String runName
        String workDir
        String revision
        String profile
        List<String> configFiles
        String paramsFile
        String entryName
        String resume
        boolean latest
        boolean stubRun
        String mainScript
        Map<String, String> params
        Launcher launcher
    }
}
