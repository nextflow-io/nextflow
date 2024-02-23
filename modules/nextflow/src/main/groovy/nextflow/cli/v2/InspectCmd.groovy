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
import groovy.util.logging.Slf4j
import nextflow.cli.CliOptions
import nextflow.cli.CmdInspect
import picocli.CommandLine.Command
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import picocli.CommandLine.ParentCommand
import picocli.CommandLine.Unmatched

/**
 * CLI `inspect` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
@Command(
    name = 'inspect',
    description = 'Inspect process settings in a pipeline project'
)
class InspectCmd extends AbstractCmd implements CmdInspect.Options {

    @ParentCommand
    private Launcher launcher

    @Parameters(description = 'Project name or repository url')
    String pipeline

    @Unmatched
    List<String> unmatched = []

    @Option(names = ['-concretize'], description = "Build the container images resolved by the inspect command")
    boolean concretize

    @Option(names = ['-c','-config'], hidden = true)
    List<String> runConfig

    @Option(names = ['-format'], description = "Inspect output format. Can be 'json' or 'config'")
    String format = 'json'

    @Option(names = ['-i','-ignore-errors'], description = 'Ignore errors while inspecting the pipeline')
    boolean ignoreErrors

    @Option(names = '-params-file', description = 'Load script parameters from a JSON/YAML file')
    String paramsFile

    @Option(names = ['-profile'], description = 'Use the given configuration profile(s)')
    String profile

    @Option(names = ['-r','-revision'], description = 'Revision of the project to inspect (either a git branch, tag or commit SHA number)')
    String revision

    Map<String,String> params

    @Override
    List<String> getArgs() { [] }

    @Override
    String getLauncherCli() {
        launcher.cliString
    }

    @Override
    CliOptions getLauncherOptions() {
        launcher.options
    }

    @Override
    void run() {
        params = ParamsHelper.parseParams(unmatched)

        final opts = new RunCmd()
        opts.launcher = launcher
        opts.ansiLog = false
        opts.preview = true
        opts.pipeline = pipeline
        opts.params = params
        opts.paramsFile = paramsFile
        opts.profile = profile
        opts.revision = revision
        opts.runConfig = runConfig

        new CmdInspect(this).run(opts)
    }

}
