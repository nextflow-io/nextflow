/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.executor.local

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.SysEnv
import nextflow.container.DockerConfig
import nextflow.executor.BashWrapperBuilder
import nextflow.processor.TaskBean
import nextflow.processor.TaskRun
import nextflow.util.CmdLineHelper
import nextflow.util.PosixPaths

/**
 * Creates the Docker command lines used by the local executor to run the
 * whole task wrapper script (i.e. {@code .command.run}) inside the task container.
 *
 * <p>Unlike {@link nextflow.container.DockerBuilder}, the commands are created as
 * argument lists that are executed without a shell on the host. This is required
 * on Windows, where the host has no Bash: the container is the only place in which
 * the task scripts can run. Host paths are mounted at the location returned by
 * {@link PosixPaths#toPosix(Path)}, the same path used in the generated scripts.
 */
@CompileStatic
class DockerWrapperLauncher {

    private final TaskBean bean

    private final DockerConfig config

    private final String name

    DockerWrapperLauncher(TaskBean bean, String name) {
        if( !bean.containerImage )
            throw new IllegalArgumentException("Missing container image for task: ${bean.name}")
        if( bean.containerConfig !instanceof DockerConfig )
            throw new IllegalArgumentException("Running the task wrapper in the container requires the Docker engine -- current engine: ${bean.containerConfig?.getEngine()}")
        this.bean = bean
        this.config = (DockerConfig) bean.containerConfig
        this.name = name
    }

    /**
     * @return The name assigned to the task container
     */
    String getName() { name }

    /**
     * @return The command to launch the task wrapper in a new container
     */
    List<String> runCommand() {
        final result = dockerCli()
        result.add('run')
        if( config.remove )
            result.add('--rm')

        // resources
        if( bean.containerCpus && config.cpuLimits )
            result.addAll(['--cpus', String.valueOf(bean.containerCpus)])
        else if( bean.containerCpus && !config.legacy )
            result.addAll(['--cpu-shares', String.valueOf(bean.containerCpus * 1024)])
        if( bean.containerCpuset )
            result.addAll(split(bean.containerCpuset))
        if( bean.containerMemory )
            result.addAll(['--memory', "${bean.containerMemory.toMega()}m".toString()])
        if( bean.containerPlatform )
            result.addAll(['--platform', bean.containerPlatform])

        // environment
        for( String it : envVars() )
            result.addAll(['-e', it])

        // mounts
        if( config.temp && config.temp != 'auto' )
            result.addAll(['-v', "${config.temp}:/tmp".toString()])
        result.addAll(volumes())
        result.addAll(['-w', PosixPaths.toPosix(bean.workDir)])

        final entrypoint = config.entrypointOverride()
        if( entrypoint )
            result.addAll(['--entrypoint', '/bin/bash'])

        // user provided options, the name comes after them so it has precedence
        result.addAll(split(config.runOptions))
        result.addAll(split(bean.containerOptions))
        result.addAll(['--name', name])
        result.add(bean.containerImage)

        // the task wrapper, executed from the task work dir
        final shell = BashWrapperBuilder.BASH
        result.addAll(entrypoint ? shell.subList(1, shell.size()) : shell)
        result.add(TaskRun.CMD_RUN)
        return result
    }

    /**
     * @return The command to stop the task container, or {@code null} when
     *      the container kill is disabled by the {@code docker.kill} setting
     */
    List<String> stopCommand() {
        final kill = config.getKill()
        if( kill == false )
            return null
        final result = dockerCli()
        if( kill instanceof CharSequence )
            result.addAll(['kill', '-s', kill.toString()])
        else
            result.add('stop')
        result.add(name)
        return result
    }

    protected List<String> dockerCli() {
        final result = new ArrayList<String>(50)
        if( config.sudo )
            result.add('sudo')
        result.add('docker')
        result.addAll(split(config.engineOptions))
        return result
    }

    protected List<String> envVars() {
        final result = new ArrayList<String>()
        // propagate the wrapper debug level -- docker copies the value from the client environment
        if( SysEnv.containsKey('NXF_DEBUG') )
            result.add('NXF_DEBUG')
        for( String it : config.getEnvWhitelist() )
            result.add(it)
        return result
    }

    /**
     * The task work dir and the directories of the input files and bin scripts
     * are mounted, avoiding nested mounts unless a read-only mount includes the
     * (always writable) task work dir.
     */
    protected List<String> volumes() {
        final readOnly = !config.writableInputMounts
        final inputs = new ArrayList<Path>()
        for( Path it : bean.inputFiles.values() ) {
            if( it.parent != null )
                inputs.add(it.parent)
        }
        if( bean.binDirs )
            inputs.addAll(bean.binDirs)
        if( bean.containerMount != null )
            inputs.add(bean.containerMount)

        final result = new ArrayList<String>()
        final workDir = bean.workDir
        if( readOnly ) {
            for( Path it : outermost(inputs) ) {
                if( !isWithin(it, workDir) )
                    result.addAll(['-v', volume(it, true)])
            }
            result.addAll(['-v', volume(workDir, false)])
        }
        else {
            inputs.add(workDir)
            for( Path it : outermost(inputs) )
                result.addAll(['-v', volume(it, false)])
        }
        return result
    }

    protected String volume(Path path, boolean readOnly) {
        final flags = new ArrayList<String>()
        if( readOnly )
            flags.add('ro')
        if( config.mountFlags )
            flags.add(config.mountFlags.trim())
        final opts = flags ? ':' + flags.join(',') : ''
        return "${path}:${PosixPaths.toPosix(path)}${opts}".toString()
    }

    /**
     * @return The given paths without duplicates and without the paths nested within another one
     */
    static protected List<Path> outermost(List<Path> paths) {
        final sorted = paths.unique(false).sort(false) { Path it -> PosixPaths.toPosix(it) }
        final result = new ArrayList<Path>()
        for( Path it : sorted ) {
            if( !result.any { Path parent -> isWithin(it, parent) } )
                result.add(it)
        }
        return result
    }

    static protected boolean isWithin(Path path, Path parent) {
        final p = PosixPaths.toPosix(path)
        final q = PosixPaths.toPosix(parent)
        return p == q || p.startsWith(q.endsWith('/') ? q : q + '/')
    }

    static protected List<String> split(String options) {
        return options ? CmdLineHelper.splitter(options) : List.<String>of()
    }

}
