/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.util

import java.util.concurrent.CompletableFuture

import nextflow.SysEnv
import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Timeout(5)
class ThreadsTest extends Specification {

    def 'should create and start platform thread' () {
        given:
        def result = new CompletableFuture()

        when:
        def thread = Threads.start { result.complete('done') }
        then:
        !thread.isVirtual()
        and:
        result.get() == 'done'
    }

    def 'should create and start platform thread with name' () {
        given:
        def result = new CompletableFuture()

        when:
        def thread = Threads.start('abc') { result.complete('done') }
        then:
        thread.name == 'abc'
        !thread.isVirtual()
        and:
        result.get() == 'done'
    }

    def 'should create and start virtual thread' () {
        given:
        SysEnv.push(NXF_ENABLE_VIRTUAL_THREADS: 'true')
        def result = new CompletableFuture()

        when:
        def thread = Threads.start{ result.complete('done') }
        then:
        thread.isVirtual()
        and:
        result.get() == 'done'

        cleanup:
        SysEnv.pop()
    }

    def 'should create and start virtual thread with name' () {
        given:
        SysEnv.push(NXF_ENABLE_VIRTUAL_THREADS: 'true')
        def result = new CompletableFuture()

        when:
        def thread = Threads.start('abc') { result.complete('done') }
        then:
        thread.name == 'abc'
        thread.isVirtual()
        and:
        result.get() == 'done'

        cleanup:
        SysEnv.pop()
    }
}
