/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.processor

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BatchContextTest extends Specification {

    def 'should add entries to collector' () {

        given:
        def c = new BatchContext()

        when:
        c.collect('job-1')
        c.collect('job-2')
        c.collect('job-2')
        then:
        c.size() == 2

    }


    def 'should return a batch of ids' () {

        given:
        def c = new BatchContext()
        (0..55).each { c.collect( "id-$it".toString() ) }

        expect:
        c.getBatchFor('id-0', 10) == 'id-0'..'id-9'
        c.getBatchFor('id-10', 10) == 'id-10'..'id-19'
        c.getBatchFor('id-50', 10) == 'id-50'..'id-55'
    }


    def 'should add some entries in the cache' (){

        given:
        def c = new BatchContext()

        when:
        c.put( 'id-1', 'Foo' )
        c.put( 'id-2', 'Bar' )
        c.put( 'id-3', 'Hello')

        then:
        c.contains('id-1')
        c.contains('id-2')
        c.contains('id-3')
        !c.contains('id-4')
        c.get('id-1') == 'Foo'
        c.get('id-3') == 'Hello'
    }
}
