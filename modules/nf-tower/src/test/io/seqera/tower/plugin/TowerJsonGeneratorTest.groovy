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

package io.seqera.tower.plugin

import groovy.json.JsonGenerator
import groovy.json.JsonSlurper
import nextflow.trace.ProgressRecord
import nextflow.trace.WorkflowStats
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TowerJsonGeneratorTest extends Specification {

    def 'should chomp too long values' () {
        given:
        def scheme = [foo: 5, 'bar.one': 5]
        def gen = new TowerJsonGenerator(new JsonGenerator.Options(), scheme)

        when:
        def x = gen.toJson( [foo: "Hola", bar: 'mundo'] )
        then:
        x == '{"foo":"Hola","bar":"mundo"}'

        when:
        x = gen.toJson( [foo: "Hello world"] )
        then:
        x == '{"foo":"Hello"}'

        when:
        x = gen.toJson( [bar: [one: "Hello world", two: "Hola mundo"]] )
        then:
        x == '{"bar":{"one":"Hello","two":"Hola mundo"}}'
    }

    def 'should normalise gitmodules attribute' () {
        given:
        def scheme = ['workflow.manifest.gitmodules': 10]
        def gen = new TowerJsonGenerator(new JsonGenerator.Options(), scheme)

        when:
        def json = gen.toJson( [workflow: [manifest: [gitmodules: ['a','b','c']]]] )
        then:
        json == '{"workflow":{"manifest":{"gitmodules":"a,b,c"}}}'

        when:
        json = gen.toJson( [workflow: [manifest: [gitmodules: 'abc']]] )
        then:
        json == '{"workflow":{"manifest":{"gitmodules":"abc"}}}'

        when:
        json = gen.toJson( [workflow: [manifest: [gitmodules: '123456789012345']]] )
        then:
        json == '{"workflow":{"manifest":{"gitmodules":"1234567890"}}}'
    }

    def 'should serialise progress records' () {
        given:
        def gen = new TowerJsonGenerator(new JsonGenerator.Options(), [:])
        and:
        def rec1 = new ProgressRecord(1, 'foo')

        and:
        def rec2 = new ProgressRecord(2, 'bar')
        rec2.pending = 1
        rec2.submitted = 2
        rec2.running = 3
        rec2.succeeded = 4
        rec2.failed = 5
        rec2.aborted = 6
        rec2.stored = 7
        rec2.ignored = 8
        rec2.retries = 9
        rec2.cached = 10
        rec2.loadCpus = 11
        rec2.loadMemory = 12
        rec2.peakRunning = 13
        rec2.peakCpus = 14
        rec2.peakMemory = 15
        rec2.terminated = true

        when:
        def json = gen.toJson([progress: [rec1, rec2]])
        then:
        def copy = (Map)new JsonSlurper().parseText(json)
        copy.size() == 1

        and:
        def progress = (List<Map>)copy.progress
        progress.size() == 2

        and:
        progress.get(0) == [
                index:1,
                name: 'foo',
                pending: 0,
                submitted: 0,
                running: 0,
                succeeded: 0,
                failed: 0,
                aborted: 0,
                stored: 0,
                ignored: 0,
                retries: 0,
                cached: 0,
                loadCpus: 0,
                loadMemory: 0,
                peakCpus: 0,
                peakMemory: 0,
                peakRunning: 0,
                terminated: false]
        and:
        progress[1] == [
                index:2,
                name: 'bar',
                pending: 1,
                submitted: 2,
                running: 3,
                succeeded: 4,
                failed: 5,
                aborted: 6,
                stored: 7,
                ignored: 8,
                retries: 9,
                cached: 10,
                loadCpus: 11,
                loadMemory: 12,
                peakRunning: 13,
                peakCpus: 14,
                peakMemory: 15,
                terminated: true]

    }


    def 'should serialise workflow stats' () {
        given:
        def gen = new TowerJsonGenerator(new JsonGenerator.Options(), [:])
        and:
        def rec1 = new ProgressRecord(1, 'foo')
        def rec2 = new ProgressRecord(2, 'bar')
        rec2.pending = 1
        rec2.submitted = 2
        rec2.running = 3
        rec2.succeeded = 4
        rec2.failed = 5
        rec2.aborted = 6
        rec2.stored = 7
        rec2.ignored = 8
        rec2.retries = 9
        rec2.cached = 10
        rec2.loadCpus = 11
        rec2.loadMemory = 12
        rec2.peakRunning = 13
        rec2.peakCpus = 14
        rec2.peakMemory = 15
        rec2.terminated = true
        and:
        def stats = new WorkflowStats(
                succeededCount: 1,
                cachedCount: 2,
                failedCount: 3,
                ignoredCount: 4,
                pendingCount: 5,
                submittedCount: 6,
                runningCount: 7,
                retriesCount: 8,
                abortedCount: 9,
                records: [1:rec1, 2:rec2])

        when:
        def json = gen.toJson(new WorkflowProgress(stats))
        then:
        def copy = (Map)new JsonSlurper().parseText(json)
        copy.succeeded == 1
        copy.cached == 2
        copy.failed == 3
        copy.ignored == 4
        copy.pending == 5
        copy.submitted == 6
        copy.running == 7
        copy.retries == 8
        copy.aborted == 9
        and:
        (copy.processes as List).size() == 2
        and:
        with(copy.processes[0] as Map) {
            index == 1
            name == 'foo'
            pending == 0
            submitted == 0
            running == 0
            succeeded == 0
            failed == 0
            aborted == 0
            stored == 0
            ignored == 0
            retries == 0
            cached == 0
            loadCpus == 0
            loadMemory == 0
            peakCpus == 0
            peakMemory == 0
            peakRunning == 0
            terminated == false
        }
        and:
        with(copy.processes[1] as Map) {
            index == 2
            name == 'bar'
            pending ==1
            submitted == 2
            running == 3
            succeeded == 4
            failed == 5
            aborted == 6
            stored == 7
            ignored == 8
            retries == 9
            cached == 10
            loadCpus == 11
            loadMemory == 12
            peakRunning == 13
            peakCpus == 14
            peakMemory == 15
            terminated == true
        }
    }

}
