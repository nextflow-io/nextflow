/*
 * Copyright 2023, Seqera Labs
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
import nextflow.cli.ILauncherOptions
import picocli.CommandLine.Option

/**
 * Top-level CLI (v2) options
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class LauncherOptions implements ILauncherOptions {

    boolean ansiLogCli

    void setAnsiLog(boolean value) { ansiLogCli = value }

    @Option(names = ['--bg'], arity = '0', description = 'Execute nextflow in background')
    boolean background

    @Option(names = ['-C'], description = 'Use the specified configuration file(s), ignoring any defaults')
    List<String> config

    @Option(names = ['-c','--config'], description = 'Add the specified file to configuration set')
    List<String> userConfig

    @Option(names = ['--config-ignore-includes'], description = 'Disable the parsing of config includes')
    boolean ignoreConfigIncludes

    @Option(names = ['-D'], description = 'Set JVM properties')
    Map<String,String> jvmOpts

    @Option(names = ['--debug'], description = 'Enable DEBUG level logging for the specified package name -- multiple packages can be provided as a comma-separated list (e.g. \'-debug nextflow,io.seqera\')', hidden = true)
    List<String> debug

    @Option(names = ['-d','--dockerize'], arity = '0', description = 'Launch Nextflow via Docker (experimental)')
    boolean dockerize

    @Option(names = ['-h'], description = 'Print this help', usageHelp = true)
    boolean help

    @Option(names = ['--log'], description = 'Set the log file path')
    String logFile

    @Option(names = ['-q','--quiet'], description = 'Do not print information messages')
    boolean quiet

    @Option(names = ['--self-update'], arity = '0', description = 'Update Nextflow to the latest version', hidden = true)
    boolean selfUpdate

    @Option(names = ['--syslog'], arity = '0..1', fallbackValue = 'localhost', description = 'Send logs to syslog server (e.g. localhost:514)')
    String syslog

    @Option(names = ['--trace'], description = 'Enable TRACE level logging for the specified package name -- multiple packages can be provided as a comma-separated list (e.g. \'-trace nextflow,io.seqera\')')
    List<String> trace

    @Option(names = ['-v'], description = 'Print the version number and exit')
    boolean version

    @Option(names = ['-V','--version'], description = 'Print the full version info and exit')
    boolean fullVersion

}
