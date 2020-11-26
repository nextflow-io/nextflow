/*
 * Copyright 2020, Seqera Labs
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
import spock.lang.Specification
import spock.lang.Unroll

import java.nio.file.Files

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdRunTest extends Specification {

    @Unroll
    def 'should parse cmd param=#STR' () {

        expect:
        CmdRun.parseParamValue(STR)  == EXPECTED

        where:
        STR         | EXPECTED
        null        | null
        'true'      | true
        'false'     | false
        'foo'       | 'foo'
        '10'        | 10i
        '3000000000'| 3000000000l
        '20.33'     | 20.33d
        '--foo'     | '--foo'
    }

    def 'should parse nested params' () {
        when:
        CmdRun.addParam(PARAMS, KEY, VALUE)
        then:
        PARAMS == EXPECTED

        where:
        PARAMS          | KEY       | VALUE     | EXPECTED
        [:]             | 'foo'     | '1'       | [foo: 1]
        [foo: 1]        | 'bar'     | '2'       | [foo: 1, bar: 2]
        [:]             | 'x.y.z'   | 'Hola'    | [x: [y: [z: 'Hola']]]
        [a: [p:1], x:3] | 'a.q'     | '2'       | [a: [p:1, q: 2], x:3]
        [:]             | /x\.y\.z/ | 'Hola'    | ['x.y.z': 'Hola']
        [:]             | /x.y\.z/  | 'Hola'    | ['x': ['y.z': 'Hola']]
    }

    def 'should return parsed config' () {
        given:
        def cmd = new CmdRun(profile: 'first', withTower: 'http://foo.com', launcher: new Launcher())
        def base = Files.createTempDirectory('test')
        base.resolve('nextflow.config').text = '''
        profiles {
            first {
                params {
                  foo = 'Hello world'
                  awsKey = 'xyz'
                }
                process {
                    executor = { 'local' }
                }
            }
            second {
                params.none = 'Blah'
            }
        }
        '''
        when:
        def txt = cmd.resolveConfig(base)
        then:
        txt == '''\
            params {
               foo = 'Hello world'
               awsKey = '[secret]'
            }
            
            process {
               executor = { 'local' }
            }

            workDir = 'work'
            
            tower {
               enabled = true
               endpoint = 'http://foo.com'
            }
            '''.stripIndent()

        cleanup:
        base?.deleteDir()
    }

    @Unroll
    def 'should check run name #STR' () {
        expect:
        CmdRun.matchRunName(STR) == EXPECTED
        where:
        EXPECTED    | STR
        true        | 'foo'
        true        | 'f00'
        true        | 'f-00'
        true        | 'f-a-b'
        true        | 'f-0-1'
        true        | 'foo-bar'
        true        | 'a' * 80
        and:
        true        | 'f_00'
        true        | 'f_a_b'
        true        | 'f_0_1'
        true        | 'foo_bar'
        and:
        false       | '0foo'
        false       | '-foo'
        false       | 'foo--bar'
        false       | 'foo__bar'
        false       | 'foo_-bar'
        false       | 'a' * 81

    }


    def 'should parse params file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def JSON = '{"abc": 1, "xyz": 2}'
        def YAML = '''
                    ---
                    foo: 1
                    bar: 2
                    '''.stripIndent()

        when:
        def file = folder.resolve('params.json')
        file.text = JSON
        and:
        def cmd = new CmdRun(paramsFile: file.toString())
        def params = cmd.getParsedParams()
        then:
        params.abc == 1
        params.xyz == 2
        
        when:
        file = folder.resolve('params.yaml')
        file.text = YAML
        and:
        cmd = new CmdRun(paramsFile: file.toString())
        params = cmd.getParsedParams()
        then:
        params.foo == 1
        params.bar == 2

        when:
        cmd = new CmdRun(env: [NXF_PARAMS_FILE: file.toString()])
        params = cmd.getParsedParams()
        then:
        params.foo == 1
        params.bar == 2


        when:
        cmd = new CmdRun(env: [NXF_PARAMS_FILE: '/missing/path'])
        cmd.getParsedParams()
        then:
        def e = thrown(AbortOperationException)
        e.message == 'Specified params file does not exists: /missing/path'


        cleanup:
        folder?.delete()
    }
}
