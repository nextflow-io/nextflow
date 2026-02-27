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

package nextflow.config.control

import java.nio.file.Path

import org.codehaus.groovy.syntax.SyntaxException
import spock.lang.Specification
import test.TestUtils

import static test.TestUtils.deleteDir
import static test.TestUtils.tempDir
import static test.TestUtils.tempFile

/**
 * @see nextflow.config.control.ResolveIncludeVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ResolveIncludeTest extends Specification {

    List<SyntaxException> check(List<Path> files) {
        final parser = new ConfigParser()
        return TestUtils.check(parser, files)
    }

    def 'should report an error for an invalid include source' () {
        given:
        def root = tempDir()
        def config = tempFile(root, 'nextflow.config',
            '''\
            includeConfig 'extra.config'
            ''')
        def extra = root.resolve('extra.config')

        when:
        def errors = check([config])
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 13
        errors[0].getOriginalMessage() == "Invalid include source: '${extra}'"

        cleanup:
        deleteDir(root)
    }

    def 'should resolve an include' () {
        given:
        def root = tempDir()
        def config = tempFile(root, 'nextflow.config',
            '''\
            includeConfig 'extra.config'
            ''')
        def extra = tempFile(root, 'extra.config',
            '''\
            process.cpus = 8
            '''
        )

        when:
        def errors = check([config, extra])
        then:
        errors.size() == 0

        cleanup:
        deleteDir(root)
    }

    def 'should not resolve remote includes' () {
        given:
        def root = tempDir()
        def config = tempFile(root, 'nextflow.config',
            '''\
            includeConfig 'http://example.com/nextflow.config'
            ''')

        when:
        def errors = check([config])
        then:
        errors.size() == 0

        cleanup:
        deleteDir(root)
    }

    def 'should not resolve dynamic includes' () {
        given:
        def root = tempDir()
        def config = tempFile(root, 'nextflow.config',
            '''\
            includeConfig ({ 'extra.config' }())
            ''')

        when:
        def errors = check([config])
        then:
        errors.size() == 0

        cleanup:
        deleteDir(root)
    }

}
