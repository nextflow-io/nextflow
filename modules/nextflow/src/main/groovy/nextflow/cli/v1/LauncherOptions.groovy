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
import groovy.transform.CompileStatic
import nextflow.cli.ILauncherOptions

/**
 * Top-level CLI (v1) options
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class LauncherOptions implements ILauncherOptions {

    boolean ansiLogCli

    void setAnsiLog(boolean value) { ansiLogCli = value }

    @Parameter(names = ['-bg'], arity = 0, description = 'Execute nextflow in background')
    boolean background

    @Parameter(names = ['-C'], description = 'Use the specified configuration file(s) overriding any defaults')
    List<String> config

    @Parameter(names = ['-c','-config'], description = 'Add the specified file to configuration set')
    List<String> userConfig

    @Parameter(names = ['-config-ignore-includes'], description = 'Disable the parsing of config includes')
    boolean ignoreConfigIncludes

    @DynamicParameter(names = ['-D'], description = 'Set JVM properties' )
    Map<String,String> jvmOpts = [:]

    @Parameter(names = ['-debug'], hidden = true)
    List<String> debug

    @Parameter(names = ['-d','-dockerize'], arity = 0, description = 'Launch nextflow via Docker (experimental)')
    boolean dockerize

    @Parameter(names = ['-h'], description = 'Print this help', help = true)
    boolean help

    @Parameter(names = ['-log'], description = 'Set nextflow log file path')
    String logFile

    @Parameter(names = ['-q','-quiet'], description = 'Do not print information messages' )
    boolean quiet

    @Parameter(names = ['-self-update'], arity = 0, description = 'Update nextflow to the latest version', hidden = true)
    boolean selfUpdate

    @Parameter(names = ['-syslog'], description = 'Send logs to syslog server (eg. localhost:514)' )
    String syslog

    @Parameter(names = ['-trace'], description = 'Enable trace level logging for the specified package name - multiple packages can be provided separating them with a comma e.g. \'-trace nextflow,io.seqera\'')
    List<String> trace

    @Parameter(names = ['-v'], description = 'Print the program version')
    boolean version

    @Parameter(names = ['-version'], description = 'Print the program version (full)')
    boolean fullVersion

}
