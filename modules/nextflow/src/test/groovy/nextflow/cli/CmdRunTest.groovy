/*
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
        CmdRun.parseParam(STR)  == EXPECTED

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

}
