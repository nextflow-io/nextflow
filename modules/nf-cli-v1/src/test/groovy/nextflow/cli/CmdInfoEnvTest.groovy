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

package nextflow.cli

import nextflow.SysEnv
import spock.lang.Specification

/**
 * Check that `nextflow info` never prints the values of secret environment variables
 */
class CmdInfoEnvTest extends Specification {

    def cleanup() {
        SysEnv.pop()
    }

    def 'should not print secret environment variables' () {
        given:
        SysEnv.push(NXF_VER: '26.10.0', NXF_GITHUB_ACCESS_TOKEN: 'nxf-canary-1', AWS_SECRET_ACCESS_KEY: 'nxf-canary-2', FOO_PASSWORD: 'nxf-canary-3', FOO: 'bar')

        when:
        def info = CmdInfo.getInfo(1)
        then:
        info.contains('NXF_VER=26.10.0')
        info.contains('NXF_GITHUB_ACCESS_TOKEN=[secret]')
        !info.contains('AWS_SECRET_ACCESS_KEY')
        !info.contains('FOO=bar')
        !info.contains('nxf-canary')

        when:
        info = CmdInfo.getInfo(2)
        then:
        info.contains('NXF_GITHUB_ACCESS_TOKEN=[secret]')
        info.contains('AWS_SECRET_ACCESS_KEY=[secret]')
        info.contains('FOO_PASSWORD=[secret]')
        info.contains('FOO=bar')
        !info.contains('nxf-canary')
    }
}
