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

package nextflow.agent

import nextflow.config.ConfigValidator
import nextflow.util.Duration
import org.junit.Rule
import spock.lang.Specification
import test.OutputCapture

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AgentConfigTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

    def 'should build from a config map'() {
        when:
        def config = new AgentConfig([defaultModel: 'openai/gpt-5-mini', maxIterationsDefault: 7, requestTimeout: '90s'])
        then:
        config.defaultModel == 'openai/gpt-5-mini'
        config.maxIterationsDefault == 7
        config.requestTimeout == Duration.of('90s')
    }

    def 'should build from an empty config map'() {
        when:
        def config = new AgentConfig([:])
        then:
        config.defaultModel == null
        config.maxIterationsDefault == null
        config.requestTimeout == null
    }

    def 'should be a recognized config scope'() {
        when:
        new ConfigValidator().validate([
            agent: [
                defaultModel: 'openai/gpt-5-mini',
                maxIterationsDefault: 7,
                requestTimeout: '90s'
            ]
        ])
        then:
        !capture.toString().contains("Unrecognized config option 'agent")
    }
}
