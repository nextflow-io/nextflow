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
import com.beust.jcommander.IStringConverter
import com.beust.jcommander.Parameter
import com.beust.jcommander.Parameters
import groovy.transform.CompileStatic
import nextflow.cli.ILauncherOptions
import nextflow.cli.RunImpl
import nextflow.util.Duration

/**
 * CLI `run` sub-command (v1)
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Parameters(commandDescription = 'Execute a pipeline project')
class RunCmd extends AbstractCmd implements RunImpl.Options, HubOptions {

    static class DurationConverter implements IStringConverter<Long> {
        @Override
        Long convert(String value) {
            if( !value ) throw new IllegalArgumentException()
            if( value.isLong() ) {  return value.toLong() }
            return Duration.of(value).toMillis()
        }
    }

    @Parameter(description = 'Project name or repository url')
    List<String> args = []

    @Parameter(names = ['-ansi'], hidden = true, arity = 0)
    void setAnsi(boolean value) {
        launcher.options.ansiLog = value
    }

    @Parameter(names = ['-ansi-log'], arity = 1, description = 'Enable/disable ANSI console logging')
    void setAnsiLog(boolean value) {
        launcher.options.ansiLog = value
    }

    @Parameter(names = ['-bg'], arity = 0, hidden = true)
    void setBackground(boolean value) {
        launcher.options.background = value
    }

    @Parameter(names = ['-bucket-dir'], description = 'Remote bucket where intermediate result files are stored')
    String bucketDir

    @Parameter(names = ['-cache'], arity = 1, description = 'Enable/disable processes caching')
    Boolean cacheable

    @DynamicParameter(names = ['-cluster.'], description = 'Set cluster options', hidden = true )
    Map<String,String> clusterOptions = [:]

    @Parameter(names = ['-c','-config'], hidden = true )
    List<String> runConfig

    @Parameter(names = ['-disable-jobs-cancellation'], description = 'Prevent the cancellation of child jobs on execution termination')
    Boolean disableJobsCancellation

    @Parameter(names = ['-dsl1'], description = 'Execute the workflow using DSL1 syntax')
    boolean dsl1

    @Parameter(names = ['-dsl2'], description = 'Execute the workflow using DSL2 syntax')
    boolean dsl2

    @Parameter(names = ['-dump-channels'], description = 'Dump channels for debugging purpose')
    String dumpChannels

    @Parameter(names = ['-dump-hashes'], description = 'Dump task hash keys for debugging purpose')
    boolean dumpHashes

    @Parameter(names = ['-entry'], arity = 1, description = 'Entry workflow name to be executed')
    String entryName

    @DynamicParameter(names = ['-e.','-env.'], description = 'Add the specified variable to execution environment')
    Map<String,String> env = [:]

    @DynamicParameter(names = ['-executor.'], description = 'Set executor options', hidden = true )
    Map<String,String> executorOptions = [:]

    @Parameter(names = ['-E'], description = 'Exports all current system environment')
    boolean exportSysEnv

    @Parameter(names = ['-latest'], description = 'Pull latest changes before run')
    boolean latest

    @Parameter(names = ['-lib'], description = 'Library extension path')
    String libPath

    @Parameter(names = ['-main-script'], description = 'The script file to be executed when launching a project directory or repository' )
    String mainScript

    @Parameter(names = ['-name'], description = 'Assign a mnemonic name to the a pipeline run')
    String runName

    @Parameter(names = ['-offline'], description = 'Do not check for remote project updates')
    boolean offline

    @Parameter(names = ['-params-file'], description = 'Load script parameters from a JSON/YAML file')
    String paramsFile

    @Parameter(names = ['-plugins'], description = 'Specify the plugins to be applied for this run e.g. nf-amazon,nf-tower')
    String plugins

    @Parameter(names = ['-pi','-poll-interval'], description = 'Executor poll interval (duration string ending with ms|s|m)', converter = DurationConverter, hidden = true)
    long pollInterval

    @Parameter(names = ['-ps','-pool-size'], description = 'Number of threads in the execution pool', hidden = true)
    Integer poolSize

    @Parameter(names = ['-preview'], description = "Run the workflow script skipping the execution of all processes")
    boolean preview

    @DynamicParameter(names = ['-process.'], description = 'Set process options' )
    Map<String,String> processOptions = [:]

    @Parameter(names = ['-profile'], description = 'Choose a configuration profile')
    String profile

    @Parameter(names = ['-qs','-queue-size'], description = 'Max number of processes that can be executed in parallel by each executor')
    Integer queueSize

    @Parameter(names = ['-resume'], description = 'Execute the script using the cached results, useful to continue executions that was stopped by an error')
    String resume

    @Parameter(names = ['-r','-revision'], description = 'Revision of the project to run (either a git branch, tag or commit SHA number)')
    String revision

    @Parameter(names = ['-stdin'], hidden = true)
    boolean stdin

    @Parameter(names = ['-stub-run','-stub'], description = 'Execute the workflow replacing process scripts with command stubs')
    boolean stubRun

    @Parameter(names = ['-test'], description = 'Test a script function with the name specified')
    String test

    @Parameter(names = ['-with-apptainer'], description = 'Enable process execution in a Apptainer container')
    String withApptainer

    @Parameter(names = ['-with-charliecloud'], description = 'Enable process execution in a Charliecloud container runtime')
    String withCharliecloud

    @Parameter(names = ['-with-conda'], description = 'Use the specified Conda environment package or file (must end with .yml|.yaml suffix)')
    String withConda

    @Parameter(names = ['-without-conda'], description = 'Disable the use of Conda environments')
    Boolean withoutConda

    @Parameter(names = ['-with-dag'], description = 'Create pipeline DAG file')
    String withDag

    @Parameter(names = ['-with-docker'], description = 'Enable process execution in a Docker container')
    String withDocker

    @Parameter(names = ['-without-docker'], arity = 0, description = 'Disable process execution with Docker')
    boolean withoutDocker

    @Parameter(names = ['-with-fusion'], hidden = true)
    String withFusion

    @Parameter(names = ['-with-mpi'], hidden = true)
    boolean withMpi

    @Parameter(names = ['-N','-with-notification'], description = 'Send a notification email on workflow completion to the specified recipients')
    String withNotification

    @Parameter(names = ['-with-podman'], description = 'Enable process execution in a Podman container')
    String withPodman

    @Parameter(names = ['-without-podman'], description = 'Disable process execution in a Podman container')
    boolean withoutPodman

    @Parameter(names = ['-with-report'], description = 'Create processes execution html report')
    String withReport

    @Parameter(names = ['-with-singularity'], description = 'Enable process execution in a Singularity container')
    String withSingularity

    @Parameter(names = ['-with-spack'], description = 'Use the specified Spack environment package or file (must end with .yaml suffix)')
    String withSpack

    @Parameter(names = ['-without-spack'], description = 'Disable the use of Spack environments')
    Boolean withoutSpack

    @Parameter(names = ['-with-timeline'], description = 'Create processes execution timeline file')
    String withTimeline

    @Parameter(names = ['-with-tower'], description = 'Monitor workflow execution with Seqera Tower service')
    String withTower

    @Parameter(names = ['-with-trace'], description = 'Create processes execution tracing file')
    String withTrace

    @Parameter(names = ['-with-wave'], hidden = true)
    String withWave

    @Parameter(names = ['-with-weblog'], description = 'Send workflow status messages via HTTP to target URL')
    String withWebLog

    @Parameter(names = ['-w', '-work-dir'], description = 'Directory where intermediate result files are stored')
    String workDir

    @DynamicParameter(names = ['--'], description = 'Pipeline parameters', hidden = true)
    Map<String,String> params = [:]

    @Override
    String getPipeline() {
        stdin ? '-' : args[0]
    }

    @Override
    List<String> getArgs() {
        args.size() > 1 ? args[1..-1] : []
    }

    @Override
    String getLauncherCliString() {
        launcher.cliString
    }

    @Override
    ILauncherOptions getLauncherOptions() {
        launcher.options
    }

    @Override
    String getName() { 'run' }

    @Override
    void run() {
        new RunImpl(this).run()
    }

}
