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
import picocli.CommandLine.Option

/**
 * Top-level CLI (v2) options
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class CliOptionsV2 extends CliOptions {

    Boolean ansiLogCli

    void setAnsiLog(boolean value) { ansiLogCli = value }

    @Option(names = ['--bg'], arity = '0', description = 'Execute nextflow in background')
    boolean background

    @Option(names = ['-C'], split = ',', description = 'Use the specified configuration file(s), overriding any defaults')
    List<String> config

    @Option(names = ['-c','--config'], split = ',', paramLabel = '<config>', description = 'Add the specified file to configuration set')
    List<String> userConfig

    @Option(names = ['--config-ignore-includes'], description = 'Disable the parsing of config includes')
    boolean ignoreConfigIncludes

    @Option(names = ['-D'], paramLabel = '<name>=<value>', description = 'Set JVM properties')
    Map<String,String> jvmOpts = [:]

    @Option(names = ['--debug'], split = ',', paramLabel = '<package>', description = 'Enable DEBUG level logging for the specified package name', hidden = true)
    List<String> debug

    @Option(names = ['--log'], paramLabel = '<file>', description = 'Set the log file path')
    String logFile

    @Option(names = ['-q','--quiet'], description = 'Do not print information messages')
    boolean quiet

    @Option(names = ['-remote-debug'], description = "Enable JVM interactive remote debugging (experimental)")
    boolean remoteDebug

    @Option(names = ['--self-update'], arity = '0', description = 'Update Nextflow to the latest version', hidden = true)
    boolean selfUpdate

    @Option(names = ['--syslog'], arity = '0..1', fallbackValue = 'localhost', paramLabel = '<url>', description = 'Send logs to syslog server (e.g. localhost:514)')
    String syslog

    @Option(names = ['--trace'], split = ',', paramLabel = '<package>', description = 'Enable TRACE level logging for the specified package name')
    List<String> trace

    @Option(names = ['-v'], description = 'Print the version number and exit')
    boolean version

    @Option(names = ['-V','--version'], description = 'Print the full version info and exit')
    boolean fullVersion

}
