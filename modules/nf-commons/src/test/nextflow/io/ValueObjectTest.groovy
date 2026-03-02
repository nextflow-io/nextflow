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

package nextflow.io

import groovy.test.GroovyAssert
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ValueObjectTest extends Specification {

    void 'should validate ValueObjectBase annotation'() {
        expect:
        GroovyAssert.assertScript '''
            import nextflow.io.SerializableMarker
            import nextflow.io.SerializableObject

            @SerializableObject
            class Foo {
                String foo
            }

            def x = new Foo(foo: 'hello')
            assert x instanceof Serializable
            '''
    }


    void 'should validate ValueObject annotation'() {
        expect:
        GroovyAssert.assertScript '''
            import nextflow.io.SerializableMarker
            import nextflow.io.ValueObject

            @ValueObject
            class Foo {
                String foo
            }

            def x = new Foo(foo: 'hello')
            def y = new Foo(foo: 'hello')
            assert x instanceof Serializable
            assert x == y
            '''
    }

    void 'should be cloneable'() {
        expect:
        GroovyAssert.assertScript '''
            import nextflow.io.SerializableMarker
            import nextflow.io.ValueObject

            @ValueObject
            class Obj {
                String foo
                String bar
            }

            def obj = new Obj(foo: 'hello', bar:'world')
            assert obj == obj.clone()
            '''
    }

    void 'should be copyable'() {
        expect:
        GroovyAssert.assertScript '''
            import nextflow.io.SerializableMarker
            import nextflow.io.ValueObject

            @ValueObject
            class Obj {
                String foo
                String bar
            }

            def obj = new Obj(foo: 'hello', bar:'world')
            def that = obj.copyWith(foo:'hola')

            assert that.foo == 'hola'
            assert that.bar == 'world'
            '''
    }

}
