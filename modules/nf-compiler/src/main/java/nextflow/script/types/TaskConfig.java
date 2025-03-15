/*
 * Copyright 2024-2025, Seqera Labs
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
package nextflow.script.types;

import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import nextflow.script.dsl.Constant;
import nextflow.script.dsl.Description;

public interface TaskConfig {

    /* PROPERTIES */

    @Constant("attempt")
    @Description("""
        The current task attempt.
    """)
    int getAttempt();

    @Constant("exitStatus")
    @Description("""
        The exit code of the task script. Only available after the task has been executed.
    """)
    int getExitStatus();

    @Constant("hash")
    @Description("""
        The unique task hash.
    """)
    String getHash();

    @Constant("id")
    @Description("""
        The pipeline-level task index.
    """)
    int getId();

    @Constant("index")
    @Description("""
        The process-level task index.
    """)
    int getIndex();

    @Constant("name")
    @Description("""
        The current task name.
    """)
    String getName();

    @Constant("previousException")
    @Description("""
        The exception reported by the previous task attempt.
    """)
    Exception getPreviousException();

    @Constant("previousTrace")
    @Description("""
        The trace record associated with the previous task attempt.
    """)
    Map<String,?> getPreviousTrace();

    @Constant("process")
    @Description("""
        The current process name.
    """)
    String getProcess();

    @Constant("workDir")
    @Description("""
        The task unique directory.
    """)
    Path getWorkDir();

    /* DIRECTIVES */

    @Constant("accelerator")
    @Description("""
        The `accelerator` directive allows you to request hardware accelerators (e.g. GPUs) for the task execution.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#accelerator)
    """)
    Map<String,?> getAccelerator();

    @Constant("afterScript")
    @Description("""
        The `afterScript` directive allows you to execute a custom (Bash) snippet *after* the task script.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#afterscript)
    """)
    String getAfterScript();

    @Constant("arch")
    @Description("""
        The `arch` directive allows you to define the CPU architecture to build the software used by the task.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#arch)
    """)
    String getArch();

    @Constant("array")
    @Description("""
        The `array` directive allows you to submit tasks as *job arrays* for executors that support it.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#array)
    """)
    Integer getArray();

    @Constant("beforeScript")
    @Description("""
        The `beforeScript` directive allows you to execute a custom (Bash) snippet *before* the task script.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#beforescript)
    """)
    String getBeforeScript();

    @Constant("cache")
    @Description("""
        The `cache` directive allows you to store the process results to a local cache. When the cache is enabled *and* the pipeline is launched with the `-resume` option, any task executions that are already cached will be re-used.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#cache)
    """)
    String getCache();

    @Constant("clusterOptions")
    @Description("""
        The `clusterOptions` directive allows the usage of any native configuration option accepted by your cluster submit command. You can use it to request non-standard resources or use settings that are specific to your cluster and not supported out of the box by Nextflow.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#clusteroptions)
    """)
    String getClusterOptions();

    @Constant("conda")
    @Description("""
        The `conda` directive allows for the definition of the process dependencies using the [Conda](https://conda.io) package manager.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#conda)
    """)
    String getConda();

    @Constant("container")
    @Description("""
        The `container` directive allows you to execute the process script in a container.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#container)
    """)
    String getContainer();

    @Constant("containerOptions")
    @Description("""
        The `containerOptions` directive allows you to specify any container execution option supported by the underlying container engine (ie. Docker, Singularity, etc). This can be useful to provide container settings only for a specific process.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#containeroptions)
    """)
    String getContainerOptions();

    @Constant("cpus")
    @Description("""
        The `cpus` directive allows you to define the number of (logical) CPUs required by each task.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#cpus)
    """)
    int getCpus();

    @Constant("debug")
    @Description("""
        The `debug` directive allows you to print the process standard output to Nextflow\'s standard output, i.e. the console. By default this directive is disabled.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#debug)
    """)
    boolean getDebug();

    @Constant("disk")
    @Description("""
        The `disk` directive allows you to define how much local disk storage the process is allowed to use.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#disk)
    """)
    MemoryUnit getDisk();

    @Constant("errorStrategy")
    @Description("""
        The `errorStrategy` directive allows you to define what to do when a task fails.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#errorstrategy)
    """)
    String getErrorStrategy();

    @Constant("executor")
    @Description("""
        The `executor` defines the underlying system where tasks are executed.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#executor)
    """)
    String getExecutor();

    @Constant("ext")
    @Description("""
        The `ext` is a special directive used for custom settings by some executors.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#ext)
    """)
    Map<String,?> getExt();

    @Constant("fair")
    @Description("""
        The `fair` directive, when enabled, guarantees that process outputs will be emitted in the order in which they were received.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#fair)
    """)
    boolean getFair();

    @Constant("label")
    @Description("""
        The `label` directive allows you to annotate a process with a mnemonic identifier of your choice.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#label)
    """)
    String getLabel();

    @Constant("machineType")
    @Description("""
        The `machineType` directive can be used to specify a predefined Google Compute Platform [machine type](https://cloud.google.com/compute/docs/machine-types) when using the [Google Batch](https://nextflow.io/docs/latest/google.html#cloud-batch) executor.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#machinetype)
    """)
    String getMachineType();

    @Constant("maxErrors")
    @Description("""
        The `maxErrors` directive allows you to specify the maximum number of times a process can fail when using the `retry` or `ignore` error strategy. By default this directive is disabled.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#maxerrors)
    """)
    int getMaxErrors();

    @Constant("maxForks")
    @Description("""
        The `maxForks` directive allows you to define the maximum number of tasks (per process) that can be executed in parallel.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#maxforks)
    """)
    Integer getMaxForks();

    @Constant("maxRetries")
    @Description("""
        The `maxRetries` directive allows you to define the maximum number of times a task can be retried when using the `retry` error strategy. By default only one retry is allowed.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#maxretries)
    """)
    int getMaxRetries();

    @Constant("maxSubmitAwait")
    @Description("""
        The `maxSubmitAwait` directives allows you to specify how long a task can remain in the submission queue. If a task remains in the queue beyond this time limit, it will fail.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#maxsubmitawait)
    """)
    Duration getMaxSubmitAwait();

    @Constant("memory")
    @Description("""
        The `memory` directive allows you to define how much memory is required by each task. Can be a string (e.g. `\'8 GB\'`) or a memory unit (e.g. `8.GB`).

        [Read more](https://nextflow.io/docs/latest/reference/process.html#memory)
    """)
    MemoryUnit getMemory();

    @Constant("module")
    @Description("""
        The `module` directive allows you to provide software dependencies to a process using [Environment Modules](http://modules.sourceforge.net/).

        [Read more](https://nextflow.io/docs/latest/reference/process.html#module)
    """)
    List<String> getModule();

    @Constant("penv")
    @Description("""
        The `penv` directive allows you to define the parallel environment to be used when submitting a parallel task to the [SGE](https://nextflow.io/docs/latest/executor.html#sge) resource manager.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#penv)
    """)
    String getPenv();

    @Constant("pod")
    @Description("""
        The `pod` directive allows you to define pod specific settings, such as environment variables, secrets, and config maps, when using the [Kubernetes](https://nextflow.io/docs/latest/kubernetes.html) executor.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#pod)
    """)
    List<?> getPod();

    @Constant("publishDir")
    @Description("""
        The `publishDir` directive allows you to publish the process output files to a directory.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#publishdir)
    """)
    List<?> getPublishDir();

    @Constant("queue")
    @Description("""
        The `queue` directive allows you to specify the queue to which jobs are submitted when using a grid executor.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#queue)
    """)
    String getQueue();

    @Constant("resourceLabels")
    @Description("""
        The `resourceLabels` directive allows you to specify custom name-value pairs which are applied to the compute resources used for the process execution.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#resourcelabels)
    """)
    Map<String,String> getResourceLabels();

    @Constant("resourceLimits")
    @Description("""
        The `resourceLimits` directive allows you to specify environment-specific limits for task resource requests.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#resourcelimits)
    """)
    Map<String,?> getResourceLimits();

    @Constant("scratch")
    @Description("""
        The `scratch` directive allows you to execute each task in a temporary directory that is local to the compute node.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#scratch)
    """)
    String getScratch();

    @Constant("secret")
    @Description("""
        The `secret` directive allows you to securely provide secrets to a process.

        [Read more](https://nextflow.io/docs/latest/secrets.html#process-directive)
    """)
    List<String> getSecret();

    @Constant("shell")
    @Description("""
        The `shell` directive allows you to define a custom shell command for process scripts. By default, script blocks are executed with `/bin/bash -ue`.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#shell)
    """)
    List<String> getShell();

    @Constant("spack")
    @Description("""
        The `spack` directive allows you to provide software dependencies using the [Spack](https://spack.io) package manager.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#spack)
    """)
    String getSpack();

    @Constant("stageInMode")
    @Description("""
        The `stageInMode` directive defines how input files are staged into the task work directory.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#stageinmode)
    """)
    String getStageInMode();

    @Constant("stageOutMode")
    @Description("""
        The `stageOutMode` directive defines how output files are staged out from the scratch directory to the task work directory when using the `scratch` directive.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#stageoutmode)
    """)
    String getStageOutMode();

    @Constant("storeDir")
    @Description("""
        The `storeDir` directive allows you to use an external directory as a *permanent* cache for process outputs.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#storedir)
    """)
    Path getStoreDir();

    @Constant("tag")
    @Description("""
        The `tag` directive allows you to associate each process execution with a custom label, so that it will be easier to identify in the log file or in a report.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#tag)
    """)
    String getTag();

    @Constant("time")
    @Description("""
        The `time` directive allows you to define how long a task is allowed to run.

        [Read more](https://nextflow.io/docs/latest/reference/process.html#time)
    """)
    Duration getTime();

}
