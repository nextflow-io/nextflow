/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
 * Helper methods to handle Docker containers
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DockerBuilder extends ContainerBuilder<DockerBuilder> {

    private boolean sudo

    private boolean remove = true

    private boolean userEmulation

    private String registry

    private String name

    private boolean tty

    private static final String USER_AND_HOME_EMULATION = '-u $(id -u) -e "HOME=${HOME}" -v /etc/passwd:/etc/passwd:ro -v /etc/shadow:/etc/shadow:ro -v /etc/group:/etc/group:ro -v $HOME:$HOME'

    private String removeCommand

    private String killCommand

    private kill = true

    private boolean legacy

    private String mountFlags0

    DockerBuilder( String name ) {
        this.image = name
    }

    @Override
    DockerBuilder params( Map params ) {
        if( !params ) return this

        if( params.containsKey('temp') )
            this.temp = params.temp

        if( params.containsKey('engineOptions') )
            addEngineOptions(params.engineOptions.toString())

        if( params.containsKey('runOptions') )
            addRunOptions(params.runOptions.toString())

        if ( params.containsKey('userEmulation') )
            this.userEmulation = params.userEmulation?.toString() == 'true'

        if ( params.containsKey('remove') )
            this.remove = params.remove?.toString() == 'true'

        if( params.containsKey('sudo') )
            this.sudo = params.sudo?.toString() == 'true'

        if( params.containsKey('tty') )
            this.tty = params.tty?.toString() == 'true'

        if( params.containsKey('entry') )
            this.entryPoint = params.entry

        if( params.containsKey('kill') )
            this.kill = params.kill

        if( params.containsKey('legacy') )
            this.legacy = params.legacy?.toString() == 'true'

        if( params.containsKey('readOnlyInputs') )
            this.readOnlyInputs = params.readOnlyInputs?.toString() == 'true'

        if( params.containsKey('mountFlags') )
            this.mountFlags0 = params.mountFlags

        return this
    }

    @Override
    DockerBuilder setName( String name ) {
        this.name = name
        return this
    }

    @Override
    DockerBuilder build(StringBuilder result) {
        assert image

        if( sudo )
            result << 'sudo '

        result << 'docker '

        if( engineOptions )
            result << engineOptions.join(' ') << ' '

        result << 'run -i '

        if( cpus ) {
            if( legacy )
                result << "--cpuset ${cpus} "
            else
                result << "--cpuset-cpus ${cpus} "
        }

        if( memory )
            result << "--memory ${memory} "

        if( tty )
            result << '-t '

        // add the environment
        appendEnv(result)

        if( temp )
            result << "-v $temp:/tmp "

        if( userEmulation )
            result << USER_AND_HOME_EMULATION << ' '

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

        // use an explicit 'docker rm' command since the --rm flag may fail. See https://groups.google.com/d/msg/docker-user/0Ayim0wv2Ls/tDC-tlAK03YJ
        if( remove && name ) {
            removeCommand = 'docker rm ' + name
            if( sudo ) removeCommand = 'sudo ' + removeCommand
        }

        if( kill )  {
            killCommand = 'docker kill '
            // if `kill` is a string it is interpreted as a the kill signal
            if( kill instanceof String ) killCommand += "-s $kill "
            killCommand += name
            // prefix with sudo if required
            if( sudo ) killCommand = 'sudo ' + killCommand
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
