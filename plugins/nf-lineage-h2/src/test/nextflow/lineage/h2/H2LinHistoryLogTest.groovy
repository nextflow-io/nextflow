/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.lineage.h2

import nextflow.lineage.config.LineageConfig
import spock.lang.Shared
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class H2LinHistoryLogTest extends Specification {

    @Shared
    H2LinStore store

    def setupSpec() {
        def uri = "jdbc:h2:mem:testdb;DB_CLOSE_DELAY=-1"
        def config = new LineageConfig([store:[location:uri]])
        store = new H2LinStore().open(config)
    }

    def cleanupSpec() {
        store.close()
    }

    def cleanup() {
        store.truncateAllTables()
    }

    def 'should write lid record' () {
        given:
        def log = store.getHistoryLog()
        def uuid = UUID.randomUUID()
        when:
        log.write('foo', uuid, '1234')
        then:
        noExceptionThrown()

        when:
        def rec = log.getRecord(uuid)
        then:
        rec.runName == 'foo'
        rec.sessionId == uuid
        rec.runLid == '1234'
    }

    def 'should update run lid' () {
        given:
        def log = store.getHistoryLog()
        def uuid = UUID.randomUUID()
        when:
        log.write('foo', uuid, '1234')
        then:
        noExceptionThrown()

        when:
        log.updateRunLid(uuid, '4444')
        then:
        noExceptionThrown()

        when:
        def rec = log.getRecord(uuid)
        then:
        rec.runName == 'foo'
        rec.sessionId == uuid
        rec.runLid == '4444'
    }

    def 'should update get records' () {
        given:
        def log = store.getHistoryLog()
        def uuid1 = UUID.randomUUID()
        def uuid2 = UUID.randomUUID()
        def uuid3 = UUID.randomUUID()
        when:
        log.write('foo1', uuid1, '1')
        log.write('foo2', uuid2, '2')
        log.write('foo3', uuid3, '3')
        then:
        noExceptionThrown()

        when:
        def all = log.getRecords()
        then:
        all.size()==3
        and:
        all[0].runName == 'foo1'
        all[0].sessionId == uuid1
        all[0].runLid == '1'
        and:
        all[1].runName == 'foo2'
        all[1].sessionId == uuid2
        all[1].runLid == '2'
        and:
        all[2].runName == 'foo3'
        all[2].sessionId == uuid3
        all[2].runLid == '3'
    }

}
