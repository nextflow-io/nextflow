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

package nextflow.executor

import nextflow.util.Duration
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ExecutorConfigTest extends Specification {


    def 'test get queue size'() {

        when:
        def config = new ExecutorConfig([ '$sge':[queueSize: 123] ])
        then:
        config.getQueueSize('sge', 1) == 123
        config.getQueueSize('xxx', 1) == 1
        config.getQueueSize(null, 1) == 1

        when:
        config = new ExecutorConfig([ queueSize: 321, '$sge':[queueSize:789] ])
        then:
        config.getQueueSize('sge', 2) == 789
        config.getQueueSize('xxx', 2) == 321
        config.getQueueSize(null, 2) == 321

        when:
        config = new ExecutorConfig([:])
        then:
        config.getQueueSize('sge', 1) == 1
        config.getQueueSize('xxx', 2) == 2
        config.getQueueSize(null, 3) == 3

    }

    def 'test get poll interval'() {

        when:
        def config = new ExecutorConfig(['$sge':[pollInterval: 345] ])
        then:
        config.getPollInterval('sge').toMillis() == 345
        config.getPollInterval('xxx').toMillis() == 1_000
        config.getPollInterval(null).toMillis() == 1_000
        config.getPollInterval(null, 2_000 as Duration).toMillis() == 2_000

        when:
        config = new ExecutorConfig([ pollInterval: 321, '$sge':[pollInterval:789] ])
        then:
        config.getPollInterval('sge').toMillis() == 789
        config.getPollInterval('xxx').toMillis() == 321
        config.getPollInterval(null).toMillis() == 321

    }

    def 'test get exit read timeout'() {

        setup:
        def config = new ExecutorConfig(['$sge':[exitReadTimeout: '5s'] ])

        expect:
        config.getExitReadTimeout('sge') == '5sec' as Duration
        config.getExitReadTimeout('lsf') == '270sec' as Duration

    }

    def 'test get queue stat interval'() {

        setup:
        def config = new ExecutorConfig(['$sge':[queueStatInterval: '4sec'] ])

        expect:
        config.getQueueStatInterval('sge') == '4sec' as Duration
        config.getQueueStatInterval('lsf') == '1min' as Duration

    }

    def 'test monitor dump interval'() {

        setup:
        def config = new ExecutorConfig(['$sge':[dumpInterval: '6sec'] ])

        expect:
        config.getMonitorDumpInterval('sge') == '6sec' as Duration
        config.getMonitorDumpInterval('lsf') == '5min' as Duration

    }

    def 'test get exec config prop'() {

        when:
        def config = new ExecutorConfig([cpus:123, queueSize:222, '$hazelcast': [queueSize:333] ])
        then:
        config.getExecConfigProp( 'hazelcast', 'cpus', null ) == 123
        config.getExecConfigProp( 'hazelcast', 'queueSize', null ) == 333
        config.getExecConfigProp( 'local', 'queueSize', null ) == 222
        config.getExecConfigProp( 'local', 'queueSize', 'beta') == 222
        config.getExecConfigProp( 'hazelcast', 'jobName', null ) ==  null
        config.getExecConfigProp( 'hazelcast', 'jobName', 'alpha') == 'alpha'
        config.getExecConfigProp( 'hazelcast', 'jobName', 'alpha', [NXF_EXECUTOR_JOBNAME:'hola']) == 'hola'
    }


}
