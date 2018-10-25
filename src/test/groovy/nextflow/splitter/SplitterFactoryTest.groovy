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
class SplitterFactoryTest extends Specification {

    def testCreateSplitter() {

        expect:
        SplitterFactory.create('text') instanceof TextSplitter
        SplitterFactory.create('fasta') instanceof FastaSplitter

    }

    def testArgsToOptions() {

        given:
        def closure = { -> 1 }

        expect:
        SplitterFactory.argsToOpt( [ ] as Object[] ) == [:]
        SplitterFactory.argsToOpt( [ closure ] as Object[] ) == [ each: closure ]
        SplitterFactory.argsToOpt( [ [x:1, y:2] ] as Object[] ) == [x:1, y:2]
        SplitterFactory.argsToOpt( [ [x:1, y:2], closure ] as Object[] ) == [x:1, y:2, each: closure]
    }

}
