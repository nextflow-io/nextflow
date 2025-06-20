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
 */

package nextflow.cli

import java.nio.file.Files
import java.nio.file.Path

import nextflow.SysEnv
import nextflow.exception.AbortOperationException
import nextflow.plugin.Plugins
import nextflow.secret.SecretsLoader
import org.junit.Rule
import spock.lang.Shared
import spock.lang.Specification
import test.OutputCapture
/**
 *
 * @author Jorge Aguilera <jorge.aguilera@seqera.io>
 */
class CmdSecretTest extends Specification {

    @Rule
    OutputCapture capture = new OutputCapture()

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
        SysEnv.push([NXF_SECRETS_FILE: secretFile.toString()])
        //required to run all test due collisions with others
        SecretsLoader.instance.reset()
    }

    def cleanupSpec() {
        SysEnv.pop()
        tempDir?.deleteDir()
    }

    def 'should validate #COMMAND does not accept #ARGUMENTS' () {
        when:
        new CmdSecret(args: [COMMAND] + ARGUMENTS).run()

        then:
        thrown(AbortOperationException)

        where:
        COMMAND | ARGUMENTS
        'list' | ['foo']
        'get' | ['']
        'get' | ['a','b']
        'set' | ['']
        'set' | ['a']
        'delete' | ['']
        'delete' | ['a','b']
    }


    def 'should print no secrets info' () {

        given:
        secretFile.delete()

        when:
        new CmdSecret(args: ['list']).run()
        def screen = capture
                .toString()
                .readLines()
                .join('\n')

        then:
        screen.indexOf("no secrets available") != -1
    }

    def 'should set a secret' () {

        given:
        secretFile.delete()

        when:
        new CmdSecret(args: ['set','foo','bar']).run()

        then:
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

        when:
        new CmdSecret(args: ['delete', 'foo']).run()

        then:
        secretFile.text.indexOf('"name": "foo"') == -1
    }

}
