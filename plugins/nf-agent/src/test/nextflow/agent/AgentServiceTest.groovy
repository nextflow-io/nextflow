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

import java.lang.reflect.Method

import spock.lang.Specification

class AgentServiceTest extends Specification {

    def 'should declare a single chat(String):String proxy method'() {
        when:
        Method[] methods = AgentService.getDeclaredMethods()

        then: 'exactly one declared method'
        methods.length == 1

        and: 'named chat, returning String, taking a single String argument'
        final m = methods[0]
        m.name == 'chat'
        m.returnType == String
        m.parameterTypes.toList() == [String]

        and: 'it is an interface with no annotations on the method (avoids @Agent clash)'
        AgentService.isInterface()
        m.annotations.length == 0
    }
}
