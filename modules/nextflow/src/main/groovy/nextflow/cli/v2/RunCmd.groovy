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
import nextflow.cli.RunImpl
import nextflow.util.Duration
import picocli.CommandLine.Command
import picocli.CommandLine.ITypeConverter
import picocli.CommandLine.Option
import picocli.CommandLine.Parameters
import picocli.CommandLine.ParentCommand

/**
 * CLI `run` sub-command (v2)
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Command(
    name = 'run',
    description = 'Execute a pipeline'
)
class RunCmd extends AbstractCmd implements RunImpl.Options, HubOptions {

    static class DurationConverter implements ITypeConverter<Long> {
        @Override
        Long convert(String value) {
            if( !value ) throw new IllegalArgumentException()
            if( value.isLong() ) { return value.toLong() }
            return Duration.of(value).toMillis()
        }
    }

    @ParentCommand
    private Launcher launcher

    @Parameters(index = '0', description = 'Project name or repository url')
    List<String> args

    @Option(names = ['--ansi-log'], arity = '1', description = 'Use ANSI logging')
    void setAnsiLog(boolean value) {
        launcher.ansiLog = value
    }

    @Option(names = ['--bg'], arity = '0', description = 'Run as a background process')
    void setBackground(boolean value) {
        launcher.background = value
    }

    @Option(names = ['--bucket-dir'], description = 'Remote bucket where intermediate result files are stored')
    String bucketDir

    @Option(names = ['--cache'], arity = '1', description = 'Enable/disable process caching')
    Boolean cacheable

    @Option(names = ['--cluster.'], description = 'Set cluster options', hidden = true)
    Map<String,String> clusterOptions = [:]

    @Option(names = ['-c','--config'], hidden = true)
    List<String> runConfig

    @Option(names = ['--disable-jobs-cancellation'], description = 'Do not cancel child jobs on workflow termination')
    Boolean disableJobsCancellation

    @Option(names = ['--dsl1'], description = 'Execute the workflow using DSL1 syntax (no longer supported)')
    boolean dsl1

    @Option(names = ['--dsl2'], description = 'Execute the workflow using DSL2 syntax')
    boolean dsl2

    @Option(names = ['--dump-channels'], description = 'Dump channels for debugging purposes')
    String dumpChannels

    @Option(names = ['--dump-hashes'], description = 'Dump task hash keys for debugging purposes')
    boolean dumpHashes

    @Option(names = ['--entry'], arity = '1', description = 'Entry workflow name to be executed')
    String entryName

    @Option(names = ['-e.','--env.'], description = 'Add the specified variable to execution environment')
    Map<String,String> env = [:]

    @Option(names = ['--executor.'], description = 'Set executor options', hidden = true)
    Map<String,String> executorOptions = [:]

    @Option(names = ['-E','--export-sys-env'], description = 'Export the current system environment')
    boolean exportSysEnv

    @Option(names = ['--latest'], description = 'Pull latest changes before run')
    boolean latest

    @Option(names = ['--lib'], description = 'Library extension path')
    String libPath

    @Option(names = ['--main-script'], description = 'The script file to be executed when launching a project directory or repository' )
    String mainScript

    @Option(names = ['--name'], description = 'Assign a mnemonic name to the pipeline run')
    String runName

    @Option(names = ['--offline'], description = 'Do not check for remote project updates')
    boolean offline = System.getenv('NXF_OFFLINE')=='true'

    @Option(names = ['--params-file'], description = 'Load script parameters from a JSON/YAML file')
    String paramsFile

    @Option(names = ['--plugins'], description = 'Specify the plugins to be applied for this run e.g. nf-amazon,nf-tower')
    String plugins

    @Option(names = ['--poll-interval'], description = 'Executor poll interval (duration string ending with ms|s|m)', converter = DurationConverter, hidden = true)
    long pollInterval

    @Option(names = ['--pool-size'], description = 'Number of threads in the execution pool', hidden = true)
    Integer poolSize

    @Option(names = ['--preview'], description = 'Run the workflow script skipping the execution of all processes')
    boolean preview

    @Option(names = ['--process.'], description = 'Set process options' )
    Map<String,String> processOptions = [:]

    @Option(names = ['--profile'], description = 'Use a configuration profile -- multiple profiles can be specified as a comma-separated list')
    String profile

    @Option(names = ['--queue-size'], description = 'Max number of processes that can be executed in parallel by each executor')
    Integer queueSize

    @Option(names = ['--resume'], description = 'Execute the script using the cached results, useful to continue executions that was stopped by an error')
    String resume

    @Option(names = ['-r','--revision'], description = 'Revision of the project to run (either a git branch, tag or commit SHA number)')
    String revision

    @Option(names = ['--stub-run'], description = 'Execute the workflow replacing process scripts with command stubs')
    boolean stubRun

    @Option(names = ['--test'], description = 'Test a script function with the name specified')
    String test

    @Option(names = ['-w','--work-dir'], description = 'Directory where intermediate result files are stored')
    String workDir

    @Option(names = ['--with-apptainer'], description = 'Enable process execution in a Apptainer container')
    def withApptainer

    @Option(names = ['--with-charliecloud'], description = 'Enable process execution in a Charliecloud container runtime')
    def withCharliecloud

    @Option(names = ['--with-conda'], description = 'Use the specified Conda environment package or file (must end with .yml|.yaml suffix)')
    String withConda

    @Option(names = ['--without-conda'], arity = '0', description = 'Disable the use of Conda environments')
    Boolean withoutConda

    @Option(names = ['--with-dag'], description = 'Create pipeline DAG file')
    String withDag

    @Option(names = ['--with-docker'], description = 'Enable process execution in a Docker container')
    def withDocker

    @Option(names = ['--without-docker'], arity = '0', description = 'Disable process execution with Docker')
    boolean withoutDocker

    @Option(names = ['--with-fusion'], hidden = true)
    String withFusion

    @Option(names = ['--with-mpi'], hidden = true)
    boolean withMpi

    @Option(names = ['-N','--with-notification'], description = 'Send a notification email on workflow completion to the specified recipients')
    String withNotification

    @Option(names = ['--with-podman'], description = 'Enable process execution in a Podman container')
    def withPodman

    @Option(names = ['--without-podman'], arity = '0', description = 'Disable process execution in a Podman container')
    def withoutPodman

    @Option(names = ['--with-report'], description = 'Create processes execution html report')
    String withReport

    @Option(names = ['--with-singularity'], description = 'Enable process execution in a Singularity container')
    def withSingularity

    @Option(names = ['--with-spack'], description = 'Use the specified Spack environment package or file (must end with .yaml suffix)')
    String withSpack

    @Option(names = ['--without-spack'], arity = '0', description = 'Disable the use of Spack environments')
    Boolean withoutSpack

    @Option(names = ['--with-timeline'], description = 'Create processes execution timeline file')
    String withTimeline

    @Option(names = ['--with-tower'], description = 'Monitor workflow execution with Seqera Tower service')
    String withTower

    @Option(names = ['--with-trace'], description = 'Create processes execution tracing file')
    String withTrace

    @Option(names = ['--with-wave'], hidden = true)
    String withWave

    @Option(names = ['--with-weblog'], description = 'Send workflow status messages via HTTP to target URL')
    String withWebLog

    @Parameters(description = 'Set pipeline parameters')
    List<String> params

    @Override
    Map<String,String> getParams() {
        Map<String,String> paramsMap = [:]

        for( int i = 0; i < params.size(); i += 2 ) {
            String key = params[i]
            String value = params[i + 1]
            paramsMap.put(key, value)
        }

        return paramsMap
    }

    @Override
    boolean getStdin() {
        args.size() > 0 && args[0] == '-'
    }

    @Override
    String getLauncherCliString() { getCliString() }

    @Override
    ILauncherOptions getLauncherOptions() { launcher }

    @Override
    Integer call() {
        new RunImpl(this).run()
        return 0
    }

}
