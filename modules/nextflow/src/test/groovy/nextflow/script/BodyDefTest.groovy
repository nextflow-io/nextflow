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

package nextflow.script

import spock.lang.Timeout
import test.Dsl2Spec

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class BodyDefTest extends Dsl2Spec {

    def 'should set script type properly' () {

        when:
        def body = new BodyDef({->'echo foo'}, 'echo foo', section)
        then:
        body.type == expected
        body.isShell == shell
        where:
        section     | expected              | shell
        'exec'      | ScriptType.GROOVY     | false
        'script'    | ScriptType.SCRIPTLET  | false
        'shell'     | ScriptType.SCRIPTLET  | true
        'workflow'  | ScriptType.GROOVY     | false
    }


    def 'should thrown illegal argument exception for invalid section' () {

        when:
        new BodyDef({}, 'echo foo', 'bar')
        then:
        thrown(IllegalArgumentException)

    }

    def 'should return empty set'() {
        when:
        def body = new BodyDef({->'echo foo'}, 'echo foo')
        then:
        body.getValNames() == [] as Set
    }

}
