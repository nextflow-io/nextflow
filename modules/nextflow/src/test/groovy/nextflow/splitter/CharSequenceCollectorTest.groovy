/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.splitter

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CharSequenceCollectorTest extends Specification {

    def 'test string buffer' () {

        given:
        def buffer = new CharSequenceCollector()
        when:
        buffer.add('alpha')
        buffer.add('-')
        buffer.add('gamma')
        buffer.add('-')
        buffer.add('delta')
        buffer.add(null)
        buffer.add('')

        then:
        buffer.toString() == 'alpha-gamma-delta'

    }

    def 'test is empty' () {

        when:
        def buffer = new CharSequenceCollector()
        then:
        !buffer.hasChunk()

        when:
        buffer.add('hello')
        then:
        buffer.hasChunk()

        when:
        buffer.nextChunk()
        then:
        !buffer.hasChunk()

    }

    def 'test reset' () {

        given:
        def buffer = new CharSequenceCollector()
        buffer.add('hello')

        when:
        buffer.nextChunk()
        then:
        !buffer.hasChunk()
        buffer.toString() == ''
    }

}
