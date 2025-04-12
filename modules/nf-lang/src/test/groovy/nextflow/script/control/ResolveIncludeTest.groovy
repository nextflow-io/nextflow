/*
 * Copyright 2024-2025, Seqera Labs
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

package nextflow.script.control

import java.nio.file.Path

import org.codehaus.groovy.syntax.SyntaxException
import spock.lang.Specification
import test.TestUtils

import static test.TestUtils.deleteDir
import static test.TestUtils.tempDir
import static test.TestUtils.tempFile

/**
 * @see nextflow.script.control.ResolveIncludeVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ResolveIncludeTest extends Specification {

    List<SyntaxException> check(List<Path> files) {
        final parser = new ScriptParser()
        return TestUtils.check(parser, files)
    }

    def 'should report an error for an invalid include source' () {
        given:
        def root = tempDir()
        def main = tempFile(root, 'main.nf',
            '''\
            include { hello } from './hello.nf'
            ''')
        def module = root.resolve('hello.nf')

        when:
        def errors = check([main])
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 13
        errors[0].getOriginalMessage() == "Invalid include source: '${module}'"

        cleanup:
        deleteDir(root)
    }

    def 'should report an error for an include source that could not be parsed' () {
        given:
        def root = tempDir()
        def main = tempFile(root, 'main.nf',
            '''\
            include { hello } from './modules/hello'
            ''')
        def module = tempFile(root, 'modules/hello/main.nf',
            '''\
            println 1 2 3
            ''')

        when:
        def errors = check([main, module])
        then:
        errors.size() == 2
        errors[0].getSourceLocator().endsWith('main.nf')
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 13
        errors[0].getOriginalMessage() == "Module could not be parsed: '${module}'"
        errors[1].getSourceLocator().endsWith('modules/hello/main.nf')
        errors[1].getStartLine() == 1
        errors[1].getStartColumn() == 23
        errors[1].getOriginalMessage() == "Unexpected input: '2'"

        cleanup:
        deleteDir(root)
    }

    def 'should report an error for an undefined include' () {
        given:
        def root = tempDir()
        def main = tempFile(root, 'main.nf',
            '''\
            include { hello } from './modules/hello'
            ''')
        def module = tempFile(root, 'modules/hello/main.nf',
            '''\
            workflow goodbye {
            }
            ''')

        when:
        def errors = check([main, module])
        then:
        errors.size() == 1
        errors[0].getSourceLocator().endsWith('main.nf')
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 13
        errors[0].getOriginalMessage() == "Included name 'hello' is not defined in module '${module}'"

        cleanup:
        deleteDir(root)
    }

    def 'should resolve an include' () {
        given:
        def root = tempDir()
        def main = tempFile(root, 'main.nf',
            '''\
            include { hello } from './modules/hello'
            ''')
        def module = tempFile(root, 'modules/hello/main.nf',
            '''\
            workflow hello {
            }
            ''')

        when:
        def errors = check([main, module])
        then:
        errors.size() == 0

        cleanup:
        deleteDir(root)
    }

}
