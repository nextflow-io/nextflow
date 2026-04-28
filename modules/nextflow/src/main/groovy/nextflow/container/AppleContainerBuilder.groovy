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
import groovy.util.logging.Slf4j

/**
 * Wrap a task execution inside an Apple container (apple/container) runtime.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AppleContainerBuilder extends ContainerBuilder<AppleContainerBuilder> {

    private boolean remove

    private boolean tty

    private String name

    private String capAdd

    private String removeCommand

    private String killCommand

    private kill = true

    AppleContainerBuilder(String name, AppleContainerConfig config) {
        this.image = name

        if( config.engineOptions )
            addEngineOptions(config.engineOptions)

        if( config.runOptions )
            addRunOptions(config.runOptions)

        if( config.temp )
            this.temp = config.temp

        this.remove = config.remove
        this.tty = config.tty

        if( !config.writableInputMounts )
            this.readOnlyInputs = true
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

        if( params.containsKey('capAdd') )
            this.capAdd = params.capAdd

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

        if( tty )
            result << '-t '

        if( cpus )
            result << "--cpus ${cpus} "

        if( memory )
            result << "-m ${memory} "

        if( platform )
            result << "--platform ${platform} "

        // environment variables
        appendEnv(result)

        if( temp )
            result << "-v $temp:/tmp "

        // volume mounts
        result << makeVolumes(mounts)
        result << '-w "$NXF_TASK_WORKDIR" '

        if( entryPoint )
            result << '--entrypoint ' << entryPoint << ' '

        if( runOptions )
            result << runOptions.join(' ') << ' '

        if( capAdd )
            result << '--cap-add ' << capAdd << ' '

        if( name )
            result << '--name ' << name << ' '

        // image is the final positional argument
        result << image

        runCommand = result.toString()

        if( remove && name ) {
            removeCommand = 'container rm ' + name
        }

        if( kill ) {
            killCommand = 'container stop '
            if( kill instanceof String ) killCommand = "container kill -s $kill "
            killCommand += name
        }

        return this
    }

    @Override
    String getRemoveCommand() { removeCommand }

    @Override
    String getKillCommand() { killCommand }
}
