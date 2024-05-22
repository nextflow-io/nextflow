/*
 * Copyright 2021, Sage-Bionetworks
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

package nextflow.secret

import java.nio.file.Files
import java.nio.file.StandardOpenOption

import nextflow.exception.AbortOperationException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class LocalSecretsProviderTest extends Specification {

    def 'should store secret file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def secretFile = folder.resolve('secrets.json')
        and:
        def provider = new LocalSecretsProvider(storeFile: secretFile)
        and:
        def s1 = new SecretImpl('foo', 'bar')
        def s2 = new SecretImpl('hello', 'world')

        when:
        provider.storeSecrets( [s1, s2] )
        then:
        secretFile.exists()
        secretFile.getPermissions() == 'rw-------'
        secretFile.text == '''\
                [
                  {
                    "name": "foo",
                    "value": "bar"
                  },
                  {
                    "name": "hello",
                    "value": "world"
                  }
                ]
                '''.stripIndent(true).rightTrim()

        cleanup:
        folder.deleteDir()
    }

    def 'should load secret file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def secretFile = folder.resolve('secrets.json');
        def json = '''
                [
                    {
                        "value": "bar",
                        "name": "foo"
                    },
                    {
                        "value": "world",
                        "name": "hello"
                    }
                ]
                '''.stripIndent()
        and:
        Files.write(secretFile, json.getBytes('utf-8'), StandardOpenOption.CREATE_NEW)
        and:
        secretFile.setPermissions('rw-------')

        when:
        def provider = new LocalSecretsProvider(storeFile: secretFile).load()
        then:
        provider.getSecret('foo') == new SecretImpl('foo', 'bar')
        provider.getSecret('hello') == new SecretImpl('hello','world')

        cleanup:
        folder.deleteDir()
    }

    def 'should fail with invalid permissions' () {
        given:
        def folder = Files.createTempDirectory('test')
        def secretFile = folder.resolve('secrets.json');
        def json = '''
                [
                    {
                        "value": "world",
                        "name": "hello"
                    }
                ]
                '''.stripIndent()
        and:
        Files.write(secretFile, json.getBytes('utf-8'), StandardOpenOption.CREATE_NEW)
        and:
        secretFile.setPermissions('rw-r-----')
        and:
        def provider = new LocalSecretsProvider(storeFile: secretFile)

        when:
        provider.loadSecrets()
        then:
        def e = thrown(AbortOperationException)
        e.message ==~ /Invalid permissions for secret store file.*/

        cleanup:
        folder.deleteDir()
    }


    def 'should put secrets' () {
        given:
        def folder = Files.createTempDirectory('test')
        def secretFile = folder.resolve('secrets.json');
        and:
        def FOO = new SecretImpl('foo', 'x')
        def BAR =  new SecretImpl('bar', 'y')
        and:
        def provider = new LocalSecretsProvider(storeFile: secretFile)

        when:
        provider.putSecret(FOO)
        provider.putSecret(BAR)
        then:
        provider.getSecret('foo') == FOO
        provider.getSecret('bar') == BAR
        and:
        provider.getSecret('other') == null
        and:
        provider.listSecretsNames() == ['bar', 'foo'] as Set

        cleanup:
        folder.deleteDir()
    }

    def 'should override secrets' () {
        given:
        def folder = Files.createTempDirectory('test')
        def secretFile = folder.resolve('secrets.json');
        and:
        def FOO1 = new SecretImpl('foo', 'x')
        def FOO2 = new SecretImpl('foo', 'y')
        and:
        def provider = new LocalSecretsProvider(storeFile: secretFile)

        when:
        provider.putSecret(FOO1)
        then:
        provider.getSecret('foo') == FOO1

        when:
        provider.putSecret(FOO2)
        then:
        provider.getSecret('foo') == FOO2

        cleanup:
        folder.deleteDir()
    }


    def 'should remove secrets' () {
        given:
        def folder = Files.createTempDirectory('test')
        def secretFile = folder.resolve('secrets.json');
        and:
        def FOO = new SecretImpl('foo','x')
        def BAR = new SecretImpl('bar', 'y')
        def provider = new LocalSecretsProvider(storeFile: secretFile)
        and:
        provider.putSecret(FOO)
        provider.putSecret(BAR)
        when:
        provider.removeSecret('bar') == BAR
        then:
        provider.getSecret('bar') == null
        and:
        provider.getSecret('foo') == FOO

        cleanup:
        folder.deleteDir()
    }

    def 'should create temp file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def secretFile = folder.resolve('secrets.json');
        and:
        def ALPHA = new SecretImpl('alpha','a')
        def AALPHA = new SecretImpl('aalpha', 'b')
        def DELTA = new SecretImpl('delta', 'd')
        def OMEGA = new SecretImpl('omega', 'o')

        def provider = new LocalSecretsProvider(storeFile: secretFile)
        and:
        provider.putSecret(ALPHA)
        provider.putSecret(AALPHA)
        provider.putSecret(DELTA)
        provider.putSecret(OMEGA)

        when:
        def file = provider.makeTempSecretsFile()
        then:
        file.permissions == 'rw-------'
        and:
        file.text == '''\
                     export aalpha="b"
                     export alpha="a"
                     export delta="d"
                     export omega="o"
                     '''.stripIndent()

        when:
        def env = provider.getSecretsEnv(['alpha','omega'])
        then:
        env == "source /dev/stdin <<<\"\$(cat <(grep -w -e 'alpha=.*' -e 'omega=.*' $file))\""

        when:
        def result = ['env', '-i', 'bash', '-c', "$env; env|sort"].execute().text
        then:
        result.count('alpha=a')==1
        result.count('omega=o')==1

        cleanup:
        folder?.deleteDir()
    }

    def 'should handle dollar in secrets' () {
        given:
        def folder = Files.createTempDirectory('test')
        def secretFile = folder.resolve('secrets.json');
        and:
        def DOLLAR1 = new SecretImpl('dollar', '$foo')
        def DOLLAR2 = new SecretImpl('dollar2', "\$foo")

        def provider = new LocalSecretsProvider(storeFile: secretFile)
        and:
        provider.putSecret(DOLLAR1)
        provider.putSecret(DOLLAR2)

        when:
        def file = provider.makeTempSecretsFile()
        then:
        file.permissions == 'rw-------'
        and:
        file.text == '''\
                     export dollar="\\\$foo"
                     export dollar2="\\\$foo"
                     '''.stripIndent()

        cleanup:
        folder?.deleteDir()
    }
}
