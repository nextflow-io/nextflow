/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.container

import java.nio.file.Path

import nextflow.executor.BashWrapperBuilder
import nextflow.util.Escape
import nextflow.util.MemoryUnit
import nextflow.util.PathTrie
/**
 * Base class for container wrapper builders
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
abstract class ContainerBuilder<V extends ContainerBuilder> {

    /**
     * Create a builder instance given the container engine
     */
    static ContainerBuilder create(String engine, String containerImage) {
        if( engine == 'docker' )
            return new DockerBuilder(containerImage)
        if( engine == 'podman' )
            return new PodmanBuilder(containerImage)
        if( engine == 'singularity' )
            return new SingularityBuilder(containerImage)
        if( engine == 'apptainer' )
            return new ApptainerBuilder(containerImage)
        if( engine == 'udocker' )
            return new UdockerBuilder(containerImage)
        if( engine == 'sarus' )
            return new SarusBuilder(containerImage)
        if( engine == 'shifter' )
            return new ShifterBuilder(containerImage)
        if( engine == 'charliecloud' )
            return new CharliecloudBuilder(containerImage)
        //
        throw new IllegalArgumentException("Unknown container engine: $engine")
    }

    final protected List env = []

    final protected List<Path> mounts = []

    protected List<String> runOptions = []

    protected List<String> engineOptions = []

    protected Integer cpus

    protected String cpuset

    protected String memory

    protected String temp

    protected String image

    protected Path workDir

    protected boolean readOnlyInputs

    @Deprecated
    protected String entryPoint

    protected String runCommand

    protected boolean mountWorkDir = true

    protected boolean privileged

    String getImage() { image }

    V addRunOptions(String str) {
        if( str )
            runOptions.add(str)
        return (V)this
    }

    V addEngineOptions(String str) {
        if( str )
            engineOptions.add(str)
        return (V)this
    }

    V setCpus(Integer value) {
        this.cpus = value
        return (V)this
    }

    V setCpuset(String value) {
        this.cpuset = value
        return (V)this
    }

    V setMemory( value ) {
        if( value instanceof MemoryUnit )
            this.memory = "${value.toMega()}m"

        else if( value instanceof String )
            this.memory = value

        else
            throw new IllegalArgumentException("Not a supported memory value")

        return (V)this
    }

    V setWorkDir( Path path ) {
        this.workDir = path
        return (V)this
    }

    V setName(String name) {
        return (V)this
    }

    V setTemp( String value ) {
        this.temp = value
        return (V)this
    }

    abstract V params( Map config )

    abstract V build(StringBuilder result)

    /**
     * @return The command string to run a container from Docker image
     */
    String getRunCommand() {
        if( !runCommand ) throw new IllegalStateException("Missing `runCommand` -- make sure `build` method has been invoked")
        runCommand
    }

    String getRunCommand(String launcher) {
        def run = getRunCommand()
        if( !run )
            throw new IllegalStateException("Missing `runCommand` -- make sure `build` method has been invoked")

        if( launcher ) {
            def result = run
            result += entryPoint ? " $entryPoint -c \"$launcher\"" : " $launcher"
            return result
        }
        return run + ' ' + launcher
    }

    String getKillCommand() { BashWrapperBuilder.KILL_CMD }

    String getRemoveCommand() { return null }

    @Deprecated
    final StringBuilder appendHelpers( StringBuilder wrapper ) {
        final result = getScriptHelpers()
        if( result ) wrapper.append(result)
        return wrapper
    }

    String getScriptHelpers() {
        return null
    }

    V build() {
        build(new StringBuilder())
    }

    V addEnv( entry ) {
        env.add(entry instanceof GString ? entry.toString() : entry)
        return (V)this
    }

    V addMount( Path path ) {
        if( path )
            mounts.add(path)
        return (V)this
    }

    V addMounts( List<Path> paths ) {
        if( paths ) for( Path it : paths )
            mounts.add(it)
        return (V)this
    }

    V addMountForInputs( Map<String,Path> inputFiles ) {
        mounts.addAll( inputFilesToPaths(inputFiles) )
        return (V)this
    }

    V addMountWorkDir(boolean flag) {
        this.mountWorkDir = flag
        return (V)this
    }

    static List<Path> inputFilesToPaths( Map<String,Path> inputFiles ) {

        List<Path> files = []
        inputFiles.each { name, storePath ->

            def path = storePath.getParent()
            if( path ) files << path

        }
        return files
    }


    protected CharSequence appendEnv( StringBuilder result ) {
        for( Object e : env ) {
            makeEnv(e, result) << ' '
        }
        return result
    }

    /**
     * Get the env command line option for the give environment definition
     *
     * @param env
     * @param result
     * @return
     */
    protected StringBuilder makeEnv( env, StringBuilder result = new StringBuilder() ) {

        if( env instanceof Map ) {
            short index = 0
            for( Map.Entry entry : env.entrySet() ) {
                if( index++ ) result << ' '
                result << ("-e \"${entry.key}=${entry.value}\"")
            }
        }
        else if( env instanceof String && env.contains('=') ) {
            result << '-e "' << env << '"'
        }
        else if( env instanceof String ) {
            result << "-e \"$env\""
        }
        else if( env ) {
            throw new IllegalArgumentException("Not a valid environment value: $env [${env.class.name}]")
        }

        return result
    }


    /**
     * Get the volumes command line options for the given list of input files
     *
     * @param mountPaths
     * @param binDir
     * @param result
     * @return
     */
    protected CharSequence makeVolumes(List<Path> mountPaths, StringBuilder result = new StringBuilder() ) {

        // add the work-dir to the list of container mounts
        final workDirStr = workDir?.toString()
        final allMounts = new ArrayList<Path>(mountPaths)
        if( workDir )
            allMounts << workDir

        // find the longest commons paths and mount only them
        final trie = new PathTrie()
        for( String it : allMounts ) { trie.add(it) }

        // when mounts are read-only make sure to remove the work-dir path
        final paths = trie.longest()
        if( readOnlyInputs && workDirStr && paths.contains(workDirStr) )
            paths.remove(workDirStr)

        for( String it : paths ) {
            if(!it) continue
            result << composeVolumePath(it,readOnlyInputs)
            result << ' '
        }

        // when mounts are read-only, make sure to include the work-dir as writable
        if( readOnlyInputs && workDir ) {
            result << composeVolumePath(workDirStr)
            result << ' '
        }

        // -- append by default the current path -- this is needed when `scratch` is set to true
        if( mountWorkDir ) {
            result << composeVolumePath('$NXF_TASK_WORKDIR')
            result << ' '
        }

        return result
    }

    protected String composeVolumePath( String path, boolean readOnly = false ) {
        return "-v ${escape(path)}:${escape(path)}${mountFlags(readOnly)}"
    }

    protected String escape(String path) {
        path.startsWith('$') ? "\"$path\"" : Escape.path(path)
    }

    protected String mountFlags(boolean readOnly) {
        readOnly ? ":ro" : ''
    }
}
