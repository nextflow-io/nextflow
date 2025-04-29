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
 */

package nextflow.lineage

import nextflow.Channel
import nextflow.extension.CH
import nextflow.lineage.model.Annotation
import nextflow.lineage.model.FileOutput

import java.nio.file.Path
import java.time.Instant
import java.time.OffsetDateTime

import nextflow.Session
import nextflow.lineage.config.LineageConfig
import nextflow.lineage.model.Checksum
import nextflow.lineage.model.DataPath
import nextflow.lineage.model.Parameter
import nextflow.lineage.model.Workflow
import nextflow.lineage.model.WorkflowOutput
import nextflow.lineage.model.WorkflowRun

import spock.lang.Specification
import spock.lang.TempDir

import java.time.ZoneOffset

import static nextflow.lineage.fs.LinPath.*

/**
 * Lineage channel extensions tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class LinChannelExImplTest extends Specification {

    @TempDir
    Path tempDir

    Path storeLocation
    Map configMap

    def setup() {
        storeLocation = tempDir.resolve("store")
        configMap = [linage: [enabled: true, store: [location: storeLocation.toString()]]]
    }

    def 'should get metadata'() {

        given:
        def uniqueId = UUID.randomUUID()
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript], "https://nextflow.io/nf-test/", "123456")
        def key = "testKey"
        def params = [new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")]
        def value1 = new WorkflowRun(workflow, uniqueId.toString(), "test_run", params)
        def outputs = [new Parameter("String", "output", "name")]
        def wfOutputs = new WorkflowOutput(OffsetDateTime.now(), "lid://testKey", outputs)
        def lidStore = new DefaultLinStore()
        def session = Mock(Session) {
            getConfig() >> configMap
        }
        lidStore.open(LineageConfig.create(session))
        lidStore.save(key, value1)
        lidStore.save("$key#output", wfOutputs)
        def channelLinExt = Spy(new LinChannelExImpl())

        when:
        def results = channelLinExt.viewLineage(session, 'lid://testKey')
        then:
        channelLinExt.getStore(session) >> lidStore
        and:
        results == value1


        when:
        results = channelLinExt.viewLineage(session, 'lid://testKey#output')
        then:
        channelLinExt.getStore(session) >> lidStore
        and:
        results == wfOutputs
    }

    def 'should return global query results' () {
        given:
        def uniqueId = UUID.randomUUID()
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(1234567), ZoneOffset.UTC)
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript],"https://nextflow.io/nf-test/", "123456" )
        def key = "testKey"
        def value1 = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [ new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")] )
        def key2 = "testKey2"
        def value2 = new FileOutput("/path/tp/file1", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, [new Annotation("key1","value1"), new Annotation("key2","value2")])
        def key3 = "testKey3"
        def value3 = new FileOutput("/path/tp/file2", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, [new Annotation("key2","value2"), new Annotation("key3","value3")])
        def key4 = "testKey4"
        def value4 = new FileOutput("/path/tp/file", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, [new Annotation("key4","value4"), new Annotation("key3","value3")])
        def lidStore = new DefaultLinStore()
        def session = Mock(Session) {
            getConfig() >> configMap
        }
        lidStore.open(LineageConfig.create(session))
        lidStore.save(key, value1)
        lidStore.save(key2, value2)
        lidStore.save(key3, value3)
        lidStore.save(key4, value4)
        def channelLinExt = Spy(new LinChannelExImpl())
        when:
        def results = CH.create()
        channelLinExt.queryLineage(session, results, [ "type":"FileOutput", "annotations.key":"key2", "annotations.value":"value2" ])
        then:
        channelLinExt.getStore(session) >> lidStore
        and:
        results.val == asUriString(key2)
        results.val == asUriString(key3)
        results.val == Channel.STOP
    }
}
