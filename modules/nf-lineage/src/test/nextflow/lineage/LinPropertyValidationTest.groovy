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
package nextflow.lineage

import spock.lang.Specification

/**
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class LinPropertyValidationTest extends Specification{

    def 'should throw exception when property does not exist'(){
        when:
            new LinPropertyValidator().validate(['value', 'not_existing'])
        then:
            def e = thrown(IllegalArgumentException)
            e.message.startsWith( "Property 'not_existing' doesn't exist in the lineage model")
    }

    def 'should not throw exception when property exist'(){
        when:
        new LinPropertyValidator().validate(['value', 'output'])
        then:
        noExceptionThrown()
    }
}
