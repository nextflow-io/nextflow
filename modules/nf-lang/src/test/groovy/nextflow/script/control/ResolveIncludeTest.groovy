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

package nextflow.script.control

import java.nio.file.Files
import java.nio.file.Path

import nextflow.module.spi.FallbackRemoteModuleResolver
import nextflow.module.spi.RemoteModuleResolverProvider
import org.codehaus.groovy.syntax.SyntaxException
import spock.lang.Specification
import spock.lang.TempDir
import test.TestUtils

import static test.TestUtils.tempFile

/**
 * @see nextflow.script.control.ResolveIncludeVisitor
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class ResolveIncludeTest extends Specification {

    @TempDir
    Path root

    def setup() {
        RemoteModuleResolverProvider.instance = new FallbackRemoteModuleResolver()
    }

    def cleanup() {
        RemoteModuleResolverProvider.reset()
    }

    List<SyntaxException> check(Path projectDir, List<Path> files) {
        final parser = new ScriptParser(projectDir)
        return TestUtils.check(parser, files)
    }

    def 'should report an error for an invalid include source' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            include { hello } from './hello.nf'
            ''')
        def module = root.resolve('hello.nf')

        when:
        def errors = check(root, [main])
        then:
        errors.size() == 1
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 13
        errors[0].getOriginalMessage() == "Invalid include source: '${module}'"
    }

    def 'should report an error for an include source that could not be parsed' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            include { hello } from './modules/hello'
            ''')
        def module = tempFile(root, 'modules/hello/main.nf',
            '''\
            println 1 2 3
            ''')

        when:
        def errors = check(root, [main, module])
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
    }

    def 'should report an error for an undefined include' () {
        given:
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
        def errors = check(root, [main, module])
        then:
        errors.size() == 1
        errors[0].getSourceLocator().endsWith('main.nf')
        errors[0].getStartLine() == 1
        errors[0].getStartColumn() == 13
        errors[0].getOriginalMessage() == "Included name 'hello' is not defined in module '${module}'"
    }

    def 'should report an error for an aliased type include' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            include { Thing as Alias } from './types.nf'
            ''')
        def module = tempFile(root, 'types.nf',
            '''\
            record Thing {
                id: String
            }
            ''')

        when:
        def errors = check(root, [main, module])
        then:
        errors.size() == 1
        errors[0].getSourceLocator().endsWith('main.nf')
        errors[0].getStartLine() == 1
        errors[0].getOriginalMessage() == 'Included types cannot be aliased'
    }

    def 'should require static typing in a script that includes a pipeline: #INCLUDE' () {
        given:
        def main = tempFile(root, 'main.nf',
            """\
            include { ${INCLUDE} } from './greet.nf'
            """)
        def module = tempFile(root, 'greet.nf',
            '''\
            nextflow.enable.types = true

            params {
                greeting: String = 'Hello'
            }

            workflow {
                println params.greeting
            }
            ''')

        when:
        def errors = check(root, [main, module])
        then:
        errors.size() == 1
        errors[0].getSourceLocator().endsWith('main.nf')
        errors[0].getOriginalMessage() == 'Including a pipeline requires `nextflow.enable.types = true` in the including script'

        where:
        INCLUDE << [ 'workflow as GREET', 'params as GreetParams' ]
    }

    def 'should require static typing in an included pipeline' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            nextflow.enable.types = true

            include { workflow as GREET } from './greet.nf'
            ''')
        def module = tempFile(root, 'greet.nf',
            '''\
            params {
                greeting: String = 'Hello'
            }

            workflow {
                println params.greeting
            }
            ''')

        when:
        def errors = check(root, [main, module])
        then:
        errors.size() == 1
        errors[0].getSourceLocator().endsWith('main.nf')
        errors[0].getOriginalMessage() == "An included pipeline must enable static typing -- set `nextflow.enable.types = true` in '${module}'"
    }

    def 'should resolve an include' () {
        given:
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
        def errors = check(root, [main, module])
        then:
        errors.size() == 0
    }

    def 'should not check arguments of a plugin function in a typed script' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            nextflow.enable.types = true

            include { samplesheetToList } from 'plugin/nf-schema'

            workflow {
                def rows = samplesheetToList('s.csv', 'schema.json')
                println rows
                println rows.size()
            }
            ''')

        when:
        def errors = check(root, [main])
        then:
        errors.size() == 0
    }

    // -- agent modules. An `agent` must bind through the ordinary include statement, in both the
    //    plain-file and the directory (`main.nf`) form. The directory form matters twice over
    //    because getLocalIncludeUri is duplicated verbatim in ResolveIncludeVisitor and
    //    ModuleResolver (lint/LSP vs run); only a test driving the include catches a divergence.

    private static final String AGENT_MODULE = '''\
        nextflow.enable.types = true

        agent reporter {
            model 'openai/gpt-4o'
            instruction 'You write QA reports.'

            input:
            sample: String

            output:
            stdout()

            prompt:
            """
            Report on ${sample}.
            """
        }
        '''.stripIndent()

    def 'should resolve an agent include from a module file' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            include { reporter } from './reporter.nf'

            workflow {
                reporter('s1')
            }
            '''.stripIndent())
        def module = tempFile(root, 'reporter.nf', AGENT_MODULE)

        when:
        def errors = check(root, [main, module])
        then:
        errors.size() == 0
    }

    def 'should resolve an agent include from a module directory' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            include { reporter } from './mods/reporter'

            workflow {
                reporter('s1')
            }
            '''.stripIndent())
        def module = tempFile(root, 'mods/reporter/main.nf', AGENT_MODULE)

        when:
        def errors = check(root, [main, module])
        then:
        errors.size() == 0
    }

    def 'should resolve an aliased agent include at the call site' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            include { reporter as qc } from './mods/reporter'

            workflow {
                qc('s1')
            }
            '''.stripIndent())
        def module = tempFile(root, 'mods/reporter/main.nf', AGENT_MODULE)

        when:
        def errors = check(root, [main, module])
        then:
        errors.size() == 0
    }

    def 'should resolve an agent include that is never called' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            include { reporter } from './mods/reporter'
            '''.stripIndent())
        def module = tempFile(root, 'mods/reporter/main.nf', AGENT_MODULE)

        when:
        def errors = check(root, [main, module])
        then:
        errors.size() == 0
    }

    def 'should report an error for an undefined agent include' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            include { nope } from './mods/reporter'
            '''.stripIndent())
        def module = tempFile(root, 'mods/reporter/main.nf', AGENT_MODULE)

        when:
        def errors = check(root, [main, module])
        then:
        errors.size() == 1
        errors[0].getSourceLocator().endsWith('main.nf')
        errors[0].getOriginalMessage() == "Included name 'nope' is not defined in module '${module}'"
    }

    def 'should resolve a pipeline include from the pipelines directory' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            nextflow.enable.types = true

            include { workflow as RNASEQ } from 'nf-core/rnaseq'
            ''')
        def pipeline = tempFile(root, 'pipelines/nf-core/rnaseq/main.nf',
            '''\
            nextflow.enable.types = true

            workflow {
                println 'hi'
            }
            ''')

        when:
        def errors = check(root, [main, pipeline])
        then:
        errors.isEmpty()
    }

    def 'should report an error for a missing pipeline' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            nextflow.enable.types = true

            include { workflow as RNASEQ } from 'nf-core/rnaseq'
            ''')

        when:
        def errors = check(root, [main])
        then:
        errors.size() == 1
        errors[0].getOriginalMessage().contains("Pipeline 'nf-core/rnaseq' not found in 'pipelines' directory")
    }

    def 'should attempt a failing pipeline include only once' () {
        given:
        def calls = 0
        RemoteModuleResolverProvider.instance = new FallbackRemoteModuleResolver() {
            @Override
            Path resolvePipeline(String name, Path baseDir) {
                calls++
                throw new IllegalStateException("cannot install ${name}")
            }
        }
        def main = tempFile(root, 'main.nf',
            '''\
            nextflow.enable.types = true

            include { workflow as RNASEQ } from 'nf-core/rnaseq'
            ''')

        when:
        def errors = check(root, [main])
        then:
        calls == 1
        errors.size() == 1
        errors[0].getOriginalMessage() == 'cannot install nf-core/rnaseq'
    }

    def 'should retry a failed pipeline include on a new analysis of the same compiler' () {
        given:
        def calls = 0
        RemoteModuleResolverProvider.instance = new FallbackRemoteModuleResolver() {
            @Override
            Path resolvePipeline(String name, Path baseDir) {
                calls++
                throw new IllegalStateException("cannot install ${name}")
            }
        }
        def main = tempFile(root, 'main.nf',
            '''\
            nextflow.enable.types = true

            include { workflow as RNASEQ } from 'nf-core/rnaseq'
            ''')
        def parser = new ScriptParser(root)
        parser.parse(main.toFile())

        when:
        parser.analyze()
        parser.analyze()
        then:
        calls == 2
    }

    def 'should not resolve a module include from the pipelines directory' () {
        given:
        def main = tempFile(root, 'main.nf',
            '''\
            include { FASTQC } from 'nf-core/fastqc'
            ''')
        Files.createDirectories(root.resolve('pipelines/nf-core/fastqc'))

        when:
        def errors = check(root, [main])
        then:
        errors.size() == 1
        errors[0].getOriginalMessage().contains("not found in 'modules' directory")
    }

    def 'should resolve a local and a module include by source string' () {
        given:
        def main = tempFile(root, 'main.nf', 'println 1')
        def local = tempFile(root, 'hello.nf', 'println 2')
        def module = tempFile(root, 'modules/nf-core/fastqc/main.nf', 'println 3')

        expect:
        ModuleResolver.getIncludeUri(main.toUri(), './hello.nf', root) == local.toUri()
        ModuleResolver.getIncludeUri(main.toUri(), 'nf-core/fastqc', root) == module.toRealPath().toUri()
    }
}
