/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.trace

import java.nio.file.Files

import nextflow.exception.AbortOperationException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TraceHelperTest extends Specification {

    def 'should overwrite a file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def path = folder.resolve('foo.txt')
        and:
        path.text = 'foo'

        when:
        def file = TraceHelper.newFileWriter(path, true, 'Test')
        file.write('Hola')
        file.close()
        then:
        path.text == 'Hola'

        cleanup:
        folder?.deleteDir()
    }

    def 'should not overwrite a file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def path = folder.resolve('foo.txt')
        and:
        path.text = 'foo'

        when:
        TraceHelper.newFileWriter(path, false, 'Test')
        then:
        def e = thrown(AbortOperationException)
        e.message == "Test file already exists: $path -- enable the 'test.overwrite' option in your config file to overwrite existing files"

        cleanup:
        folder?.deleteDir()
    }

}
