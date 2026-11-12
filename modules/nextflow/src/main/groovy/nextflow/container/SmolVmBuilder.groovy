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
import nextflow.util.MemoryUnit

/**
 * Wrap a task execution inside a smolvm microVM via `smolvm machine run`.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SmolVmBuilder extends ContainerBuilder<SmolVmBuilder> {

    private boolean tty

    private boolean network

    SmolVmBuilder(String name, SmolVmConfig config) {
        this.image = name

        if( config.engineOptions )
            addEngineOptions(config.engineOptions)

        if( config.runOptions )
            addRunOptions(config.runOptions)

        if( config.temp )
            this.temp = config.temp

        this.tty = config.tty
        this.network = config.network

        if( !config.writableInputMounts )
            this.readOnlyInputs = true
    }

    SmolVmBuilder(String name) {
        this(name, new SmolVmConfig([:]))
    }

    @Override
    SmolVmBuilder params( Map params ) {
        return this
    }

    @Override
    SmolVmBuilder build(StringBuilder result) {
        assert image

        result << 'smolvm '

        if( engineOptions )
            result << engineOptions.join(' ') << ' '

        result << 'machine run -i '

        if( tty )
            result << '-t '

        if( network )
            result << '--net '

        if( cpus )
            result << "--cpus ${cpus} "

        if( memory )
            result << "--mem ${MemoryUnit.of(memory).toMega()} "

        if( platform )
            result << "--oci-platform ${platform} "

        // environment variables
        appendEnv(result)

        if( temp )
            result << "-v $temp:/tmp "

        // volume mounts
        result << makeVolumes(mounts)
        result << '-w "$NXF_TASK_WORKDIR" '

        if( runOptions )
            result << runOptions.join(' ') << ' '

        result << '--image ' << image

        runCommand = result.toString()
        return this
    }
}
