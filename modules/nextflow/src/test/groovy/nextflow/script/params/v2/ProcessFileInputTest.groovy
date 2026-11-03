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

package nextflow.script.params.v2

import spock.lang.Specification

/**
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ProcessFileInputTest extends Specification {

    def 'should resolve file pattern'() {
        given:
        def filePattern
        def input

        when:
        input = new ProcessFileInput(null, null)
        then:
        input.getFilePattern([:]) == '*'

        when:
        input = new ProcessFileInput('*.txt', null)
        then:
        input.getFilePattern([:]) == '*.txt'

        when:
        filePattern = { -> "${id}.txt" }
        input = new ProcessFileInput(filePattern, null)
        then:
        input.getFilePattern([id: 'sample1']) == 'sample1.txt'
    }

    def 'should resolve file value'() {
        given:
        def value
        def input

        when:
        input = new ProcessFileInput(null, null)
        then:
        input.resolve([:]) == null

        when:
        value = { -> fastq }
        input = new ProcessFileInput(null, value)
        then:
        input.resolve([fastq: 'sample1.fastq']) == 'sample1.fastq'
    }

}
