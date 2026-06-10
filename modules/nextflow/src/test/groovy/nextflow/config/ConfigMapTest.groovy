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

package nextflow.config


import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigMapTest extends Specification {

    def 'should allow cfg attributes' () {

        given:
        def cfg = new ConfigMap([foo:1, bar:'two'])
        expect:
        cfg.foo == 1
        cfg.bar == 'two'
        cfg.unknown == null
        and:
        cfg.get('foo') == 1
        cfg.get('bar') == 'two'
        cfg.get('unknown') == null

    }

}
