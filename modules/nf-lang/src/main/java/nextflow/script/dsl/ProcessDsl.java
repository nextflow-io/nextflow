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
package nextflow.script.dsl;

import java.nio.file.Path;
import java.util.List;
import java.util.Map;

import nextflow.script.types.Duration;
import nextflow.script.types.MemoryUnit;
import nextflow.script.types.TaskConfig;

/**
 * DSL scope for process definitions.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
public interface ProcessDsl extends DslScope {

    @Constant("task")
    @Description("""
        Map of task properties, including directive values.
    """)
    TaskConfig getTask();

    @Description("""
        Define a script template.
    """)
    Path template(String path);

    interface DirectiveDsl extends DslScope {

        @Description("""
            The `accelerator` directive allows you to request hardware accelerators (e.g. GPUs) for the task execution.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#accelerator)
        """)
        void accelerator(Map<String,?> value);

        @Description("""
            The `afterScript` directive allows you to execute a custom (Bash) snippet *after* the task script.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#afterscript)
        """)
        void afterScript(String value);

        @Description("""
            The `arch` directive allows you to define the CPU architecture to build the software used by the task.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#arch)
        """)
        void arch(String value);

        @Description("""
            The `array` directive allows you to submit tasks as *job arrays* for executors that support it.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#array)
        """)
        void array(Integer value);

        @Description("""
            The `beforeScript` directive allows you to execute a custom (Bash) snippet *before* the task script.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#beforescript)
        """)
        void beforeScript(String value);

        @Description("""
            The `cache` directive allows you to store the process results to a local cache. When the cache is enabled *and* the pipeline is launched with the `-resume` option, any task executions that are already cached will be re-used.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#cache)
        """)
        void cache(String value);

        @Description("""
            The `clusterOptions` directive allows the usage of any native configuration option accepted by your cluster submit command. You can use it to request non-standard resources or use settings that are specific to your cluster and not supported out of the box by Nextflow.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#clusteroptions)
        """)
        void clusterOptions(String value);

        @Description("""
            The `conda` directive allows for the definition of the process dependencies using the [Conda](https://conda.io) package manager.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#conda)
        """)
        void conda(String value);

        @Description("""
            The `container` directive allows you to execute the process script in a container.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#container)
        """)
        void container(String value);

        @Description("""
            The `containerOptions` directive allows you to specify any container execution option supported by the underlying container engine (ie. Docker, Singularity, etc). This can be useful to provide container settings only for a specific process.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#containeroptions)
        """)
        void containerOptions(String value);

        @Description("""
            The `cpus` directive allows you to define the number of (logical) CPUs required by each task.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#cpus)
        """)
        void cpus(Integer value);

        @Description("""
            The `debug` directive allows you to print the process standard output to Nextflow\'s standard output, i.e. the console. By default this directive is disabled.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#debug)
        """)
        void debug(boolean value);

        @Description("""
            The `disk` directive allows you to define how much local disk storage the process is allowed to use.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#disk)
        """)
        void disk(MemoryUnit value);

        @Description("""
            The `errorStrategy` directive allows you to define what to do when a task fails.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#errorstrategy)
        """)
        void errorStrategy(String value);

        @Description("""
            The `executor` defines the underlying system where tasks are executed.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#executor)
        """)
        void executor(String value);

        @Description("""
            The `ext` is a special directive used for custom settings by some executors.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#ext)
        """)
        void ext(Map<String,?> value);

        @Description("""
            The `fair` directive, when enabled, guarantees that process outputs will be emitted in the order in which they were received.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#fair)
        """)
        void fair(boolean value);

        @Description("""
            The `label` directive allows you to annotate a process with a mnemonic identifier of your choice.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#label)
        """)
        void label(String value);

        @Description("""
            The `machineType` directive can be used to specify a predefined Google Compute Platform [machine type](https://cloud.google.com/compute/docs/machine-types) when using the [Google Batch](https://nextflow.io/docs/latest/google.html#cloud-batch) executor.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#machinetype)
        """)
        void machineType(String value);

        @Description("""
            The `maxErrors` directive allows you to specify the maximum number of times a process can fail when using the `retry` or `ignore` error strategy. By default this directive is disabled.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#maxerrors)
        """)
        void maxErrors(int value);

        @Description("""
            The `maxForks` directive allows you to define the maximum number of tasks (per process) that can be executed in parallel.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#maxforks)
        """)
        void maxForks(Integer value);

        @Description("""
            The `maxRetries` directive allows you to define the maximum number of times a task can be retried when using the `retry` error strategy. By default only one retry is allowed.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#maxretries)
        """)
        void maxRetries(int value);

        @Description("""
            The `maxSubmitAwait` directives allows you to specify how long a task can remain in the submission queue. If a task remains in the queue beyond this time limit, it will fail.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#maxsubmitawait)
        """)
        void maxSubmitAwait(Duration value);

        @Description("""
            The `memory` directive allows you to define how much memory is required by each task. Can be a string (e.g. `\'8 GB\'`) or a memory unit (e.g. `8.GB`).

            [Read more](https://nextflow.io/docs/latest/reference/process.html#memory)
        """)
        void memory(MemoryUnit value);

        @Description("""
            The `module` directive allows you to provide software dependencies to a process using [Environment Modules](http://modules.sourceforge.net/).

            [Read more](https://nextflow.io/docs/latest/reference/process.html#module)
        """)
        void module(String value);

        @Description("""
            The `penv` directive allows you to define the parallel environment to be used when submitting a parallel task to the [SGE](https://nextflow.io/docs/latest/executor.html#sge) resource manager.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#penv)
        """)
        void penv(String value);

        @Description("""
            The `pod` directive allows you to define pod specific settings, such as environment variables, secrets, and config maps, when using the [Kubernetes](https://nextflow.io/docs/latest/kubernetes.html) executor.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#pod)
        """)
        void pod(List<?> value);

        @Description("""
            The `publishDir` directive allows you to publish the process output files to a directory.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#publishdir)
        """)
        void publishDir(List<?> value);

        @Description("""
            The `queue` directive allows you to specify the queue to which jobs are submitted when using a grid executor.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#queue)
        """)
        void queue(String value);

        @Description("""
            The `resourceLabels` directive allows you to specify custom name-value pairs which are applied to the compute resources used for the process execution.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#resourcelabels)
        """)
        void resourceLabels(Map<String,?> value);

        @Description("""
            The `resourceLimits` directive allows you to specify environment-specific limits for task resource requests.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#resourcelimits)
        """)
        void resourceLimits(Map<String,?> value);

        @Description("""
            The `scratch` directive allows you to execute each task in a temporary directory that is local to the compute node.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#scratch)
        """)
        void scratch(String value);

        @Description("""
            The `secret` directive allows you to securely provide secrets to a process.

            [Read more](https://nextflow.io/docs/latest/secrets.html#process-directive)
        """)
        void secret(String value);

        @Description("""
            The `shell` directive allows you to define a custom shell command for process scripts. By default, script blocks are executed with `/bin/bash -ue`.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#shell)
        """)
        void shell(String value);

        @Description("""
            The `spack` directive allows you to provide software dependencies using the [Spack](https://spack.io) package manager.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#spack)
        """)
        void spack(String value);

        @Description("""
            The `stageInMode` directive defines how input files are staged into the task work directory.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#stageinmode)
        """)
        void stageInMode(String value);

        @Description("""
            The `stageOutMode` directive defines how output files are staged out from the scratch directory to the task work directory when using the `scratch` directive.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#stageoutmode)
        """)
        void stageOutMode(String value);

        @Description("""
            The `storeDir` directive allows you to use an external directory as a *permanent* cache for process outputs.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#storedir)
        """)
        void storeDir(String value);

        @Description("""
            The `tag` directive allows you to associate each process execution with a custom label, so that it will be easier to identify in the log file or in a report.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#tag)
        """)
        void tag(String value);

        @Description("""
            The `time` directive allows you to define how long a task is allowed to run.

            [Read more](https://nextflow.io/docs/latest/reference/process.html#time)
        """)
        void time(Duration value);

    }

    interface InputDsl extends DslScope {

        @Description("""
            Declare a variable input. The received value can be any type, and it will be made available to the process body (i.e. `script`, `shell`, `exec`) as a variable with the given name.
        """)
        void val(Object arg);

        @Deprecated
        @Description("""
            Declare a file input.
        """)
        void file(Object arg);

        @Description("""
            Declare a file input. The received value should be a file or collection of files.

            The argument can be an identifier or string. If an identifier, the received value will be made available to the process body as a variable. If a string, the received value will be staged into the task directory under the given alias.
        """)
        void path(Object arg);

        @Description("""
            Declare an environment variable input. The received value should be a string, and it will be exported to the task environment as an environment variable given by `identifier`.
        """)
        void env(Object arg);

        @Description("""
            Declare a `stdin` input. The received value should be a string, and it will be provided as the standard input (i.e. `stdin`) to the task script. It should be declared only once for a process.
        """)
        void stdin();

        @Description("""
            Declare a tuple input. Each argument should be an input declaration such as `val`, `path`, `env`, or `stdin`.

            The received value should be a tuple with the same number of elements as the `tuple` declaration, and each received element should be compatible with the corresponding `tuple` argument. Each tuple element is treated the same way as if it were a standalone input.
        """)
        void tuple(Object... args);

        @Deprecated
        @Description("""
            Declare an `each` input.
        """)
        void each(Object arg);

    }

    interface OutputDsl extends DslScope {

        @Description("""
            Declare a value output. The argument can be any value, and it can reference any output variables defined in the process body (i.e. variables declared without the `def` keyword).
        """)
        void val(Object arg);

        @Deprecated
        @Description("""
            Declare a file output.
        """)
        void file(Object arg);

        @Description("""
            Declare a file output. It receives the output files from the task environment that match the given pattern.
        """)
        void path(Object arg);

        @Description("""
            Declare an environment variable output. It receives the value of the environment variable given by `identifier` from the task environment.
        """)
        void env(Object arg);

        @Description("""
            Declare a `stdout` output. It receives the standard output of the task script.
        """)
        void stdout();

        @Description("""
            Declare an `eval` output. It receives the standard output of the given command, which is executed in the task environment after the task script.
        """)
        void eval(Object arg);

        @Description("""
            Declare a tuple output. Each argument should be an output declaration such as `val`, `path`, `env`, `stdin`, or `eval`. Each tuple element is treated the same way as if it were a standalone output.
        """)
        void tuple(Object... args);

    }

}
