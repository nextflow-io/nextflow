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

    def 'should remove an api key' () {
        given: 'the `api_?key` alternation of SECRET_KEYS. `accessKey` does NOT match `apiKey`, so'
        // without it a resolved LLM provider credential (`agent.apiKey`, possibly from a
        // `secrets.*` reference already expanded to its real value) is persisted VERBATIM by the
        // lineage observer and shipped to Seqera Platform as `workflow.configText`.
        // NOTE cross-cutting by design: this also hides every other `*apiKey*`/`*API_KEY*` config
        // key -- `ncbi.apiKey`, `registry.apiKey`, a user's own -- in `nextflow config` output and
        // in every lineage record. Intended: they are all credentials.
        def obj = [
                agent: [apiKey: 'sk-agent-1234', baseUrl: 'http://localhost:8000/v1'],
                ncbi: [apiKey: 'ncbi-abcd'],
                nested: [ [OPENAI_API_KEY: 'sk-env-5678'], [api_key: 'snake-9999'] ],
                keyword: 'not a credential' ]

        expect:
        SecretHelper.hideSecrets(obj) == [
                agent: [apiKey: '[secret]', baseUrl: 'http://localhost:8000/v1'],
                ncbi: [apiKey: '[secret]'],
                nested: [ [OPENAI_API_KEY: '[secret]'], [api_key: '[secret]'] ],
                keyword: 'not a credential' ]

        and: 'the endpoint is not a secret and stays readable -- it is diagnostics, not a credential'
        SecretHelper.hideSecrets([agent: [baseUrl: 'http://localhost:8000/v1']]) ==
                [agent: [baseUrl: 'http://localhost:8000/v1']]
    }

    def 'should mask an api key in an environment string too, not only in a config map' () {
        given: 'SECRET_KEYS gained `api_?key` and SECRET_REGEX had not, so the OUT-OF-BAND channel'
        // the agent docs still recommend for a containerized runner -- `env { OPENAI_API_KEY = ... }`
        // and `agent.containerOptions = '-e OPENAI_API_KEY'` -- was masked in the config dump and
        // NOT in the `NAME=value` environment lines a trace record carries
        expect:
        SecretHelper.secureEnvString('OPENAI_API_KEY=sk-1234') == 'OPENAI_API_KEY=[secure]'
        SecretHelper.secureEnvString('ANTHROPIC_API_KEY=sk-ant-1234') == 'ANTHROPIC_API_KEY=[secure]'
        SecretHelper.secureEnvString('NXF_AGENT_API_KEY=sk-nxf') == 'NXF_AGENT_API_KEY=[secure]'
        SecretHelper.secureEnvString('api_key=snake') == 'api_key=[secure]'
        SecretHelper.secureEnvString('apiKey=camel') == 'apiKey=[secure]'

        and: 'mixed with the lines that already matched, one per line'
        SecretHelper.secureEnvString('''\
                foo=hello
                OPENAI_API_KEY=sk-abcd
                git_token=909s-ds-'''
                .stripIndent() ) ==
                '''\
                foo=hello
                OPENAI_API_KEY=[secure]
                git_token=[secure]'''.stripIndent()

        and: 'a variable that merely mentions an api is untouched -- this masks credentials, not URLs'
        SecretHelper.secureEnvString('API_URL=https://api.openai.com/v1') == 'API_URL=https://api.openai.com/v1'
    }

}
