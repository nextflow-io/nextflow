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
package nextflow.container

import groovy.transform.CompileStatic

/**
 * Helper methods to handle Apple containers.
 *
 * Wraps tasks in {@code container run} commands using the Apple container CLI,
 * which runs linux/arm64 images natively on Apple Silicon via the Virtualization
 * framework. For amd64 images, {@code --rosetta} is added automatically.
 *
 * @author Joon-Klaps
 */
@CompileStatic
class AppleContainerBuilder extends ContainerBuilder<AppleContainerBuilder> {

    private boolean remove

    private String registry

    private String name

    private String removeCommand

    private String killCommand

    private Object kill = true

    AppleContainerBuilder(String name, AppleContainerConfig config) {
        this.image = name

        if( config.engineOptions )
            addEngineOptions(config.engineOptions)

        this.remove = config.remove

        if( config.registry )
            this.registry = config.registry

        if( config.runOptions )
            addRunOptions(config.runOptions)

        if( config.temp )
            this.temp = config.temp
    }

    AppleContainerBuilder(String name) {
        this(name, new AppleContainerConfig([:]))
    }

    @Override
    AppleContainerBuilder params( Map params ) {
        if( !params ) return this

        if( params.containsKey('entry') )
            this.entryPoint = params.entry

        if( params.containsKey('kill') )
            this.kill = params.kill

        if( params.containsKey('readOnlyInputs') )
            this.readOnlyInputs = params.readOnlyInputs?.toString() == 'true'

        return this
    }

    @Override
    AppleContainerBuilder setName( String name ) {
        this.name = name
        return this
    }

    @Override
    AppleContainerBuilder build(StringBuilder result) {
        assert image

        result << 'container '

        if( engineOptions )
            result << engineOptions.join(' ') << ' '

        result << 'run -i '

        // add the environment
        appendEnv(result)

        if( temp )
            result << "-v $temp:/tmp "

        // mount the input folders
        result << makeVolumes(mounts)
        result << '-w "$NXF_TASK_WORKDIR" '

        if( entryPoint )
            result << '--entrypoint ' << entryPoint << ' '

        if( runOptions )
            result << runOptions.join(' ') << ' '

        if( cpus )
            result << "--cpus ${cpus} "

        if( memory )
            result << "--memory ${memory.toUpperCase()} "

        // --platform selects the image variant; --rosetta enables x86_64 execution via Rosetta.
        // Both are needed for amd64 images: without --platform the CLI defaults to linux/arm64.
        if( platform ) {
            result << "--platform ${platform} "
            if( !platform.contains("arm64") )
                result << '--rosetta '
        }

        // the name is after the user option so it has precedence over any options provided by the user
        if( name )
            result << '--name ' << name << ' '

        if( registry )
            result << registry

        // finally the container image
        result << image

        // return the run command as result
        runCommand = result.toString()

        // use an explicit 'container rm' command after the container stops
        if( remove && name ) {
            removeCommand = 'container rm ' + name
        }

        if( kill && name ) {
            killCommand = 'container stop '
            // if `kill` is a string it is interpreted as the kill signal
            if( kill instanceof String ) killCommand = "container kill --signal $kill "
            killCommand += name
        }

        return this
    }

    @Override
    String getRunCommand(String launcher) {
        if( !launcher )
            return getRunCommand()
        def result = getRunCommand()
        result += entryPoint ? " -c \"$launcher\"" : " $launcher"
        return result
    }

    /**
     * @return The command string to remove a container
     */
    @Override
    String getRemoveCommand() { removeCommand }

    /**
     * @return The command string to stop/kill a running container
     */
    @Override
    String getKillCommand() { killCommand }

}
