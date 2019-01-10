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

package nextflow.extension

import spock.lang.Specification

import nextflow.Channel

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CountLinesOpTest extends Specification {

    def 'should count lines' () {

        when:
        def str = '''
            line 1
            line 2
            line 3
            line 4
            line 5
            '''
                .stripIndent().strip()

        def str2 = '''
            line 6
            line 7
            line 8
            '''
                .stripIndent().strip()

        def str3 = '''
            line 9
            line 10
            line 11
            '''
                .stripIndent().strip()

        def result = Channel.from(str, str2, str3).countLines()
        then:
        result.val == 11

    }
}
