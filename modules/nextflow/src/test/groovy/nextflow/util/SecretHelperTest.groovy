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

package nextflow.util


import test.BaseSpec
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SecretHelperTest extends BaseSpec {

    def 'should remove secret key' () {
        expect:
        SecretHelper.secureEnvString('a=b') == 'a=b'
        and:
        SecretHelper.secureEnvString('aws_key=12345') == 'aws_key=[secure]'
        SecretHelper.secureEnvString('AWS_KEY=12345') == 'AWS_KEY=[secure]'

        SecretHelper.secureEnvString('''\
                foo=hello
                aws_key=d7sds89
                git_token=909s-ds-'''
                .stripIndent() ) ==
                '''\
                foo=hello
                aws_key=[secure]
                git_token=[secure]'''.stripIndent()

    }

    def 'should remove secrets' () {
        given:
        def obj = [foo: 'hello',
                   awsKey: '1234',
                   githubToken: 'abc' ]

        expect:
        SecretHelper.hideSecrets(obj) == [
                foo: 'hello',
                awsKey: '[secret]',
                githubToken: '[secret]' ]
    }

    def 'should remove nested secrets' () {
        given:
        def obj = [
                aws: [secretKey: 'abc', accessKey: 'zzz'],
                github: [ [token: 'xxx'], [endpoint: 'this is good'] ]
                ]

        expect:
        SecretHelper.hideSecrets(obj) == [
                aws: [accessKey: '[secret]', secretKey: '[secret]'],
                github: [ [token: '[secret]'], [endpoint: 'this is good']
                ]
        ]
    }

}
