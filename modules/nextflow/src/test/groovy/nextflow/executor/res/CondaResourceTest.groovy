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

package nextflow.executor.res


import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CondaResourceTest extends Specification {

    def 'should create conda resource' () {
        expect:
        CondaResource.ofCondaPackages('foo bar').packages == 'foo bar'
        CondaResource.ofPipPackages('foo bar').pip == 'foo bar'
        CondaResource.of([packages: 'foo', pip: 'bar']).packages == 'foo'
        CondaResource.of([packages: 'foo', pip: 'bar']).pip == 'bar'
    }

    def 'should validate equals and hashcode' () {
        given:
        def c1 = CondaResource.of([packages: 'foo', pip: 'bar'])
        def c2 = CondaResource.of([packages: 'foo', pip: 'bar'])
        def c3 = CondaResource.of([packages: 'foo', pip: 'baz'])

        expect:
        c1 == c2
        c1 != c3
        and:
        c1.hashCode() == c2.hashCode()
        c1.hashCode() != c3.hashCode()
    }

    def 'should return extended syntax' () {
        expect:
        CondaResource.of(VALUE).getExtendedPackages() == EXPECTED

        where:
        VALUE                           | EXPECTED
        [packages:'foo']                | 'foo'
        [packages:'foo bar']            | 'foo bar'
        and:
        [pip:'foo']                     | 'pip:foo'
        [pip:'foo bar']                 | 'pip:foo pip:bar'
        and:
        [packages:'one', pip:'two three' ]    | 'one pip:two pip:three'

    }

    def "should check if it's a file" () {
        expect:
        CondaResource.isFilePath(VALUE)  == EXPECTED
        where:
        VALUE                   | EXPECTED
        null                    | false
        '/foo'                  | true
        'this/that'             | true
        'foo.yaml'              | true
        'foo.yml'               | true
        'http://foo.com'        | true
        'https://foo.com'       | true
        and:
        'foo bar'               | false
    }

    def 'should report and error' () {
        when:
        CondaResource.of([:])
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Conda directive cannot be empty'

        when:
        CondaResource.of([foo:'value'])
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Not a valid Conda directive attribute: 'foo'"

        when:
        CondaResource.of([packages:'/file.txt', pip:'x y'])
        then:
        e = thrown(IllegalArgumentException)
        e.message == "Conda environment file and Pip packages conflict each other - offending enviroment: '/file.txt'"
    }
}
