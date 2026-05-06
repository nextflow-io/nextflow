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

package nextflow.plugin.extension

import java.nio.file.Path

import nextflow.Channel
import nextflow.exception.DuplicateModuleFunctionException
import nextflow.exception.MissingProcessException
import nextflow.plugin.Plugins
import nextflow.plugin.TestPluginManager
import spock.lang.Shared
import spock.lang.TempDir
import test.Dsl2Spec

import static test.ScriptHelper.*
/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class PluginExtensionMethodsTest extends Dsl2Spec {

    @TempDir
    @Shared
    Path folder

    @Shared String pluginsMode

    def setup() {
        // reset previous instances
        PluginExtensionProvider.reset()
        // this need to be set *before* the plugin manager class is created
        pluginsMode = System.getProperty('pf4j.mode')
        System.setProperty('pf4j.mode', 'dev')
        // the plugin root should
        def root = Path.of('.').toAbsolutePath().normalize()
        def manager = new TestPluginManager(root)
        Plugins.init(root, 'dev', manager)
    }

    def cleanup() {
        Plugins.stop()
        PluginExtensionProvider.reset()
        pluginsMode ? System.setProperty('pf4j.mode',pluginsMode) : System.clearProperty('pf4j.mode')
    }

    def 'should execute custom operator extension/1' () {
        given:
        def SCRIPT_TEXT = '''
            include { goodbye } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.of('Bye bye folks').goodbye()
            }
            '''

        when:
        def result = runScript(SCRIPT_TEXT)

        then:
        result.val == 'Bye bye folks'
        result.val == Channel.STOP
    }

    def 'should execute custom operator extension/2' () {
        given:
        def SCRIPT_TEXT = '''
            include { reverse; goodbye } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.of('Bye bye folks').goodbye()
            }
            '''
        when:
        def result = runScript(SCRIPT_TEXT)

        then:
        result.val == 'Bye bye folks'
        result.val == Channel.STOP
    }

    def 'should execute custom factory extension/1' () {
        given:
        def SCRIPT_TEXT = '''
            include { reverse } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.reverse('a string')
            }
            '''

        when:
        def result = runScript(SCRIPT_TEXT)

        then:
        result
        result.val == 'a string'.reverse()
        result.val == Channel.STOP

    }

    def 'should execute custom factory extension/2' () {
        given:
        def SCRIPT_TEXT = '''
            include { reverse } from 'plugin/nf-test-plugin-hello'
            include { goodbye } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.reverse('a string')
            }
            '''

        when:
        def result = runScript(SCRIPT_TEXT)

        then:
        result
        result.val == 'a string'.reverse()
        result.val == Channel.STOP

    }

    def 'should execute custom operator as alias extension' () {
        given:
        def SCRIPT_TEXT = '''
            include { goodbye as myFunction } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.of(100,200,300).myFunction()
            }
            '''

        when:
        def result = runScript(SCRIPT_TEXT)

        then:
        result.val == 100
        result.val == 200
        result.val == 300
        result.val == Channel.STOP
    }

    def 'should execute custom factory as alias extension' () {
        given:
        def SCRIPT_TEXT = '''
            include { reverse as myFunction } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.myFunction('reverse this string')
            }
            '''

        when:
        def result = runScript(SCRIPT_TEXT)

        then:
        result
        result.val == 'reverse this string'.reverse()
        result.val == Channel.STOP

    }

    def 'should not include operators without the right signature' () {
        given:
        def SCRIPT_TEXT = '''
            include { goodbyeWrongSignature } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.of('Bye bye folks') | goodbyeWrongSignature
            }
            '''

        when:
        runScript(SCRIPT_TEXT)

        then:
        thrown(MissingProcessException)

    }

    def 'should not include factories without the right signature' () {
        given:
        def SCRIPT_TEXT = '''
            include { reverseCantBeImportedBecauseWrongSignature } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.reverseCantBeImportedBecauseWrongSignature('a string')
            }
            '''

        when:
        runScript(SCRIPT_TEXT)

        then:
        thrown(IllegalStateException)

    }

    def 'should execute custom functions'() {
        when:
        def result = runScript(SCRIPT_TEXT)

        then:
        result.val == EXPECTED
        result.val == Channel.STOP

        where:
        SCRIPT_TEXT                                                                           | EXPECTED
        "include { sayHello } from 'plugin/nf-test-plugin-hello'; workflow { channel.of( sayHello() ) }"     | 'hi'
        "include { sayHello } from 'plugin/nf-test-plugin-hello'; workflow { channel.of( sayHello('es') ) }" | 'hola'
        "include { sayHello as hi } from 'plugin/nf-test-plugin-hello'; workflow { channel.of( hi() ) }"     | 'hi'

    }

    def 'should call init plugin in custom functions'() {
        when:
        def result = runScript("""
            include { sayHello } from 'plugin/nf-test-plugin-hello'

            workflow {
                sayHello()
            }
            """)

        then:
        noExceptionThrown()
    }

    def 'should throw function not found'() {
        given:
        def SCRIPT_TEXT = '''
            include { sayHelloNotExist } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.of( sayHelloNotExist() )
            }
            '''

        when:
        runScript(SCRIPT_TEXT)

        then:
        thrown(IllegalStateException)
    }

    def 'should not allow to include an existing function'() {
        given:
        def SCRIPT_TEXT = '''
            def sayHello() { 'hi' }

            include { sayHello } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.of( sayHello() )
            }
            '''

        when:
        runScript(SCRIPT_TEXT)

        then:
        def e = thrown(Exception)
        e.cause.message.contains('`sayHello` is already included')
    }

    def 'should allows to include an existing function but as alias'() {
        given:
        def SCRIPT_TEXT= '''
            def sayHello() { 'hi' }

            include { sayHello as anotherHello } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.of( anotherHello() )
            }
            '''

        when:
        def result = runScript(SCRIPT_TEXT)

        then:
        result.val == 'hi'
        result.val == Channel.STOP
    }

    def 'should not include a non annotated function'() {
        given:
        def SCRIPT_TEXT= '''
            include { aNonImportedFunction } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.of( aNonImportedFunction() )
            }
            '''

        when:
        runScript(SCRIPT_TEXT)

        then:
        thrown(IllegalStateException)
    }

    def 'should execute a function in a module'() {
        given:
        def SCRIPT = folder.resolve('main.nf')
        def MODULE = folder.resolve('module.nf')

        MODULE.text = '''
        include { sayHello } from 'plugin/nf-test-plugin-hello'

        process foo {
            input:
              val lng
            output:
              stdout

            script:
            "${sayHello(lng)}"
        }
        '''

        SCRIPT.text = '''
        include { foo } from './module.nf'

        workflow {
            foo( 'en' )
        }
        '''

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'hi'

    }

    def 'should execute a function in two modules'() {
        given:
        def SCRIPT = folder.resolve('main.nf')
        def MODULE1 = folder.resolve('module1.nf')
        def MODULE2 = folder.resolve('module2.nf')

        MODULE1.text = '''
        include { sayHello } from 'plugin/nf-test-plugin-hello'

        process foo {
            input:
              val lng
            output:
              stdout

            script:
            "${sayHello(lng)}"
        }
        '''

        MODULE2.text = '''
        include { sayHello } from 'plugin/nf-test-plugin-hello'

        process bar {
            input:
              val lng
            output:
              stdout

            script:
            "${sayHello('es')}"
        }
        '''

        SCRIPT.text = '''
        include { foo } from './module1.nf'
        include { bar } from './module2.nf'

        workflow {
            foo( 'en' ) | bar
        }
        '''

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'hola'
    }

    def 'should execute custom functions and channel extension at the same time'() {
        given:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = '''
            include { sayHello; goodbye } from 'plugin/nf-test-plugin-hello'

            workflow {
                channel.of( sayHello() ).goodbye()
            }
            '''

        when:
        def result = runScript(SCRIPT)

        then:
        result.val == 'hi'
        result.val == Channel.STOP
    }

    def 'should not allow a function with the same name as a process'() {
        given:
        def SCRIPT = folder.resolve('main.nf')

        SCRIPT.text = """
        include { sayHello } from 'plugin/nf-test-plugin-hello'

        process sayHello {
            input:
              val lng
            output:
              stdout

            script:
            "Hi"
        }

        workflow {
            channel.of('hi') | sayHello
        }
        """.stripIndent()

        when:
        runScript(SCRIPT)

        then:
        def e = thrown(Exception)
        e.cause.message.contains('`sayHello` is already included')
    }

    def 'should not allow a function and a process with the same from other module'() {
        given:
        def SCRIPT = folder.resolve('main.nf')
        def MODULE1 = folder.resolve('module1.nf')

        MODULE1.text = '''
        process sayHello {
            input:
              val lng
            output:
              stdout

            script:
            "$lng"
        }
        '''
        SCRIPT.text = """
        include { sayHello } from 'plugin/nf-test-plugin-hello'
        include { sayHello } from './module1.nf'

        process foo {
            input:
              val lng
            output:
              stdout

            script:
            "Hi"
        }

        workflow {
            channel.of('hi') | foo
        }
        """.stripIndent()

        when:
        runScript(SCRIPT)

        then:
        def e = thrown(Exception)
        e.cause.message.contains('`sayHello` is already included')
    }

}
