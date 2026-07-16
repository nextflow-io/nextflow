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

import java.lang.management.ManagementFactory

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessHelperTest extends Specification {

    def 'should get process pid' () {

        given:
        def process = new ProcessBuilder().command(['bash','-c','echo $$']).start()

        when:
        def pid = ProcessHelper.pid(process)
        and:
        process.waitFor()
        then:
        process.text.trim() == pid.toString()
    }

    def 'should return self pid' () {
        when:
        def pid = ProcessHelper.selfPid()
        then:
        pid > 0
        pid == Long.parseLong(ManagementFactory.getRuntimeMXBean().getName().split("@")[0])
    }
}
