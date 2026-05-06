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

package nextflow.plugin

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PluginRefTest extends Specification {

    def 'should parse a plugin ref' () {

        when:
        def ref = PluginRef.parse(FQID)
        then:
        ref.id == ID
        ref.version == VER


        where:
        FQID        | ID    | VER
        'foo@1.0'   | 'foo' | '1.0'
    }

    def 'should compare refs' () {
        given:
        def ref1 = new PluginRef('nf-alpha','1.0.0')
        def ref2 = new PluginRef('nf-alpha','1.0.0')
        def ref3 = new PluginRef('nf-delta','2.0.0')

        expect:
        ref1 == ref2
        ref1 != ref3
        and:
        ref1 < ref3
        and:
        ref1.hashCode() == ref2.hashCode()
        ref1.hashCode() != ref3.hashCode()

    }

}
