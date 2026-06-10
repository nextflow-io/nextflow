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
package nextflow.util

import java.nio.file.Files
import java.nio.file.Path

import nextflow.exception.AbortOperationException
import nextflow.script.dsl.Nullable
import nextflow.script.types.Bag
import nextflow.script.types.Record
import org.codehaus.groovy.runtime.typehandling.GroovyCastException
import spock.lang.Specification

/**
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypeHelperTest extends Specification {

    // helper classes

    static class Sample implements Record {
        String name
        Integer count
        @Nullable String optional
    }

    static class Params implements Record {
        List<Sample> samples
    }

    // ---- getRawType ----

    def 'should return the Class itself as raw type'() {
        expect:
        TypeHelper.getRawType(String)  == String
        TypeHelper.getRawType(Integer) == Integer
    }

    def 'should extract raw type from ParameterizedType'() {
        given:
        def type = Params.getDeclaredField('samples').getGenericType()

        expect:
        TypeHelper.getRawType(type) == List
    }

    // ---- isCollectionType ----

    def 'should determine whether a type is a collection type'() {
        expect:
        TypeHelper.isCollectionType(List)
        TypeHelper.isCollectionType(Set)
        TypeHelper.isCollectionType(Bag)
        !TypeHelper.isCollectionType(String)
        !TypeHelper.isCollectionType(Integer)
        !TypeHelper.isCollectionType(Path)
        !TypeHelper.isCollectionType(Record)
        !TypeHelper.isCollectionType(Sample)
    }

    // ---- isRecordType ----

    def 'should determine whether a type is a record type'() {
        expect:
        TypeHelper.isRecordType(Record)
        TypeHelper.isRecordType(Sample)
        !TypeHelper.isRecordType(String)
        !TypeHelper.isRecordType(Integer)
        !TypeHelper.isRecordType(Path)
        !TypeHelper.isRecordType(List)
        !TypeHelper.isRecordType(Set)
        !TypeHelper.isRecordType(Bag)
    }

    // ---- asType ----

    def 'should convert raw value to target type'() {
        given:
        def inputFile = Files.createTempFile('test', '.txt')

        expect: 'primitive type'
        TypeHelper.asType(null, String) == null
        TypeHelper.asType('42', String) == '42'
        TypeHelper.asType(42, Integer) == 42

        when: 'list type'
        def result = TypeHelper.asType([1, 2, 3], List)
        then:
        result instanceof List
        result == [1, 2, 3]

        when: 'set type'
        result = TypeHelper.asType([1, 2, 2], Set)
        then:
        result instanceof Set
        result == [1, 2] as Set

        when: 'path type'
        result = TypeHelper.asType(inputFile.toString(), Path)
        then:
        result instanceof Path
        result == inputFile

        when: 'record type'
        result = TypeHelper.asType([name: 'Alice', count: 5], Sample)
        then:
        result instanceof RecordMap
        result.name == 'Alice'
        result.count == 5

        cleanup:
        inputFile?.delete()
    }

    def 'should convert raw data structure to lists and records'() {
        when:
        def params = [
            samples: [
                [name: 'Alice', count: 3, extra: 'value']
            ]
        ]
        def result = TypeHelper.asRecordType(params, Params)
        then:
        result instanceof RecordMap
        result.samples instanceof List
        result.samples[0] instanceof RecordMap
        result.samples[0].name == 'Alice'
        result.samples[0].count == 3
        result.samples[0].optional == null
        result.samples[0].extra == 'value'
    }

    // ---- asCollectionType ----

    def 'should convert raw collection to collection type'() {
        when:
        def result = TypeHelper.asCollectionType([1, 2, 2], Bag)
        then:
        result == new HashBag<>([1, 2, 2])

        when:
        result = TypeHelper.asCollectionType([3, 1, 2], List)
        then:
        result instanceof List
        result == [3, 1, 2]

        when:
        result = TypeHelper.asCollectionType([1, 2, 2], Set)
        then:
        result instanceof Set
        result == [1, 2] as Set
    }

    def 'should convert each collection element to the element type if specified'() {
        given:
        def type = Params.getDeclaredField('samples').getGenericType()

        when:
        def result = TypeHelper.asCollectionType([ [name: 'Alice', count: 3] ], type)
        then:
        result instanceof List
        result.size() == 1
        result[0] instanceof RecordMap
        result[0].name == 'Alice'
        result[0].count == 3
    }

    // ---- asRecordType ----

    def 'should convert raw map to a record'() {
        when:
        def result = TypeHelper.asRecordType([name: 'Alice', count: 3, extra: 'value'], Sample)
        then:
        result instanceof RecordMap
        result.name == 'Alice'
        result.count == 3
        result.optional == null
        result.extra == 'value'
    }

    def 'should report error when a required field is null'() {
        when:
        TypeHelper.asRecordType([name: null, count: 5], Sample)
        then:
        thrown(AbortOperationException)
    }

    def 'should report error when a required field is absent'() {
        when:
        TypeHelper.asRecordType([name: 'Alice'], Sample)
        then:
        thrown(AbortOperationException)
    }

}
