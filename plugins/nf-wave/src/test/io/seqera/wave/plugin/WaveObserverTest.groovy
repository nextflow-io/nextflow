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
 *
 */

package io.seqera.wave.plugin

import nextflow.Session
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WaveObserverTest extends Specification {

    def 'should render containers config' () {
        given:
        def sess = Mock(Session) {getConfig() >> [:] }
        def observer = new WaveObserver(sess)
        and:
        Map<String,String> containers = [:]
        containers.foo = 'quay.io/ubuntu:latest'
        containers.bar = 'quay.io/alpine:latest'

        when:
        def result = observer.renderContainersConfig(containers)
        then:
        result == '''\
            process { withName: 'foo' { container='quay.io/ubuntu:latest' }}
            process { withName: 'bar' { container='quay.io/alpine:latest' }}
            '''.stripIndent(true)
    }

}
