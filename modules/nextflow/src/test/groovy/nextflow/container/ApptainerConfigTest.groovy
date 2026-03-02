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

import nextflow.util.Duration
import spock.lang.Specification

/**
 * Apptainer config tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class ApptainerConfigTest extends Specification {

    def 'should set default values empty map'() {
        given:
        def config = new ApptainerConfig([:])

        expect:
        config.envWhitelist == []
        config.pullTimeout.toMillis() == 20 * 60 * 1000 //20 min
    }

    def 'should create config with full map'(){
        given:
        def configMap = [
            autoMounts: false,
            cacheDir: 'cacheDir',
            enabled: true,
            engineOptions: '-q -v',
            envWhitelist: 'ENV_1,ENV_2',
            libraryDir: 'libraryDir',
            noHttps: false,
            ociAutoPull: false,
            pullTimeout: '50s',
            registry: 'http://registry.com',
            runOptions: '--contain --writable'
        ]
        def config = new ApptainerConfig(configMap)

        expect:
        config.autoMounts == false
        config.cacheDir == 'cacheDir'
        config.enabled == true
        config.engineOptions == '-q -v'
        config.envWhitelist == ['ENV_1','ENV_2']
        config.libraryDir == 'libraryDir'
        config.noHttps == false
        config.ociAutoPull == false
        config.registry == 'http://registry.com'
        config.runOptions == '--contain --writable'
        config.pullTimeout.toMillis() == 50_000 // 50s

    }
}
