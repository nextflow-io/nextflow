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

package nextflow.serde

import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneOffset

import groovy.transform.EqualsAndHashCode
import nextflow.serde.gson.GsonEncoder
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GsonEncoderTest extends Specification {

    @EqualsAndHashCode
    static class Foo {
        String name
        Instant timestamp
        OffsetDateTime datetime
    }

    def 'should serialize-deserialize an object' () {
        given:
        def encoder = new GsonEncoder<Foo>() { }
        def ts = Instant.ofEpochSecond(1742638384)
        def dt = ts.atOffset(ZoneOffset.UTC)
        def foo = new Foo(name:'Yo!', timestamp: ts, datetime: dt)
        when:
        def json = encoder.encode(foo)
        then:
        json == '{"name":"Yo!","timestamp":"2025-03-22T10:13:04Z","datetime":"2025-03-22T10:13:04Z"}'
        encoder.decode(json) == foo
    }

    def 'should encode and decode polymorphic class/1'() {
        given:
        def encoder = new MyEncoder()
        def dog = new Dog("bau", 10)
        when:
        def json = encoder.encode(dog)
        then:
        json == '{"@type":"Dog","name":"bau","barkVolume":10}'

        when:
        def animal = encoder.decode(json)
        then:
        animal == dog
    }

    def 'should encode and decode polymorphic class/1'() {
        given:
        def encoder = new MyEncoder()
        def dog = new Cat("bau", true)
        when:
        def json = encoder.encode(dog)
        then:
        json == '{"@type":"Cat","name":"bau","likesSun":true}'

        when:
        def animal = encoder.decode(json)
        then:
        animal == dog
    }
}
