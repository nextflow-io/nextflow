/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.cli

import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.secret.LocalSecretsProvider
import nextflow.secret.SecretsLoader
import spock.lang.Shared
import spock.lang.Specification

import java.nio.file.Files
import java.nio.file.Path

/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class CmdSecretTest extends Specification {

    @Shared
    Path tempDir

    @Shared
    File secretFile

    def cleanup(){
        Plugins.stop()
    }

    def setupSpec(){
        tempDir = Files.createTempDirectory('test').toAbsolutePath()
        secretFile = new File("$tempDir/store.json")
        def processEnvironmentClass = System.getenv().getClass()
        def field = processEnvironmentClass.getDeclaredField('m')
        field.accessible = true
        def map = (Map<String, String>) field.get(System.getenv())
        map.put('NXF_SECRETS_FILE', secretFile.toString())
    }

    def cleanupSpec() {
        tempDir?.deleteDir()
    }

    def 'should validate #COMMAND doesnt accept #ARGUMENTS' () {
        given:
        def buffer = new ByteArrayOutputStream()

        when:
        new CmdSecret(args: [COMMAND] + ARGUMENTS, out: new PrintStream(buffer)).run()

        then:
        thrown(AbortOperationException)

        where:
        COMMAND | ARGUMENTS
        'list' | ['foo']
        'get' | ['']
        'get' | ['a','b']
        'put' | ['']
        'put' | ['a']
        'delete' | ['']
        'delete' | ['a','b']
    }


    def 'should print no secrets info' () {

        given:
        secretFile.delete()
        def buffer = new ByteArrayOutputStream()

        when:
        new CmdSecret(args: ['list'], out: new PrintStream(buffer)).run()
        def screen = buffer.toString()

        then:
        screen == "no secrets available\n"
    }

    def 'should set a secret' () {

        given:
        secretFile.delete()
        def buffer = new ByteArrayOutputStream()

        when:
        new CmdSecret(args: ['put','foo','bar'], out: new PrintStream(buffer)).run()
        def screen = buffer.toString()

        then:
        screen == ""

        and:
        secretFile.text.indexOf('"name": "foo"') != -1
        secretFile.text.indexOf('"value": "bar"') != -1
    }

    def 'should delete a secret' () {

        given:
        secretFile.text = """[
          {
            "name": "foo",
            "value": "bar"
          }
        ]
        """
        secretFile.permissions = 'rw-------'

        def buffer = new ByteArrayOutputStream()

        when:
        new CmdSecret(args: ['delete', 'foo'], out: new PrintStream(buffer)).run()
        def screen = buffer.toString()

        then:
        screen == ""

        and:
        secretFile.text.indexOf('"name": "foo"') == -1
    }

}
