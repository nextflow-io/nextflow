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
 *
 */

package nextflow.util

import java.time.Instant

import groovy.transform.Canonical
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GsonHelperTest extends Specification{

    @Canonical
    static class Person {
        String name
        Instant birthday
    }

    def 'should ser-deserialize a custom object' () {
        given:
        def person = new Person(name: 'My name', birthday: Instant.now())
        when:
        def json = GsonHelper.toJson(person)
        then:
        GsonHelper.fromJson(json,Person) == person
    }

}
