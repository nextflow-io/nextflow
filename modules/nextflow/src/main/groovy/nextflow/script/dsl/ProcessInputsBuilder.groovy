/*
 * Copyright 2013-2023, Seqera Labs
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

package nextflow.script.dsl

import groovy.util.logging.Slf4j
import nextflow.script.params.*
import nextflow.script.ProcessConfig

/**
 * Process inputs builder DSL for the {@code ProcessFn} annotation.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
class ProcessInputsBuilder {

    private ProcessConfig config

    ProcessInputsBuilder(ProcessConfig config) {
        this.config = config
    }

    void env( Object obj ) {
        new EnvInParam(config).bind(obj)
    }

    void file( Object obj ) {
        new FileInParam(config).bind(obj)
    }

    void path( Map opts=null, Object obj ) {
        new FileInParam(config)
                .setPathQualifier(true)
                .setOptions(opts)
                .bind(obj)
    }

    void stdin( Object obj = null ) {
        def result = new StdInParam(config)
        if( obj )
            result.bind(obj)
        result
    }

    ProcessConfig getConfig() {
        return config
    }

}
