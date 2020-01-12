/*
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

package nextflow.container
import groovy.transform.CompileStatic
/**
 * Helper methods to handle podman containers
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author tbugfinder <github@online.ms>
 */
@CompileStatic
class PodmanBuilder extends ContainerBuilder<PodmanBuilder> {

    private boolean remove = true

    private String registry

    private String name

    private String removeCommand

    private String killCommand

    private kill = true

    private String mountFlags0

    PodmanBuilder( String name ) {
        this.image = name
    }

    @Override
    PodmanBuilder params( Map params ) {
        if( !params ) return this

        if( params.containsKey('temp') )
            this.temp = params.temp

        if( params.containsKey('engineOptions') )
            addEngineOptions(params.engineOptions.toString())

        if( params.containsKey('runOptions') )
            addRunOptions(params.runOptions.toString())

        if ( params.containsKey('remove') )
            this.remove = params.remove?.toString() == 'true'

        if( params.containsKey('entry') )
            this.entryPoint = params.entry

        if( params.containsKey('kill') )
            this.kill = params.kill

        if( params.containsKey('readOnlyInputs') )
            this.readOnlyInputs = params.readOnlyInputs?.toString() == 'true'

        if( params.containsKey('mountFlags') )
            this.mountFlags0 = params.mountFlags

        return this
    }

    @Override
    PodmanBuilder setName( String name ) {
        this.name = name
        return this
    }

    @Override
    PodmanBuilder build(StringBuilder result) {
        assert image

        result << 'podman '

        if( engineOptions )
            result << engineOptions.join(' ') << ' '

        result << 'run -i '

        // add the environment
        appendEnv(result)

        if( temp )
            result << "-v $temp:/tmp "

        // mount the input folders
        result << makeVolumes(mounts)
        result << ' -w "$PWD" '

        if( entryPoint )
            result << '--entrypoint ' << entryPoint << ' '

        if( runOptions )
            result << runOptions.join(' ') << ' '

        // the name is after the user option so it has precedence over any options provided by the user
        if( name )
            result << '--name ' << name << ' '

        if( registry )
            result << registry

        // finally the container name
        result << image

        // return the run command as result
        runCommand = result.toString()

        // use an explicit 'podman rm' command as --rm does not remove images in all cases.
        if( remove && name ) {
            removeCommand = 'podman rm ' + name
        }

        if( kill )  {
            killCommand = 'podman kill '
            // if `kill` is a string it is interpreted as a the kill signal
            if( kill instanceof String ) killCommand += "-s $kill "
            killCommand += name
        }

        return this
    }

    protected String mountFlags(boolean readOnly) {
        def result = super.mountFlags(readOnly)
        if( !mountFlags0 )
            return result

        result ? "${result},${mountFlags0.trim()}" : ":${mountFlags0.trim()}"
    }

    @Override
    String getRunCommand(String launcher) {
        if( launcher ) {
            def result = getRunCommand()
            result += entryPoint ? " -c \"$launcher\"" : " $launcher"
            return result
        }
        return getRunCommand() + ' ' + launcher
    }

    /**
     * @return The command string to remove a container
     */
    @Override
    String getRemoveCommand() { removeCommand }

    /**
     * @return The command string to kill a running container
     */
    @Override
    String getKillCommand() { killCommand }

}
