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

package groovy.runtime.metaclass

import spock.lang.Specification

import java.nio.file.Files
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NextflowDelegatingMetaClassTest extends Specification {


    /*
     * 'isEmpty' method is implemented by the NextflowDelegatingMetaClass class
     */
    def testEmptyAndIsEmpty() {

        when:
        def file = File.createTempFile('hello','file')
        then:
        file.empty()
        file.isEmpty()
        when:
        file.text = 'hello'
        then:
        !file.empty()
        !file.isEmpty()

        when:
        def path = Files.createTempFile('hello','path')
        then:
        path.empty()
        path.isEmpty()
        when:
        path.text = 'hello'
        then:
        !path.empty()
        !path.isEmpty()

        cleanup:
        file?.delete()
        path?.delete()
    }



}
