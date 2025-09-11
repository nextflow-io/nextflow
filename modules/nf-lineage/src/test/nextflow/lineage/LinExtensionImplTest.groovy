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

import java.nio.file.Path
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneOffset

import nextflow.Channel
import nextflow.Session
import nextflow.extension.CH
import nextflow.lineage.config.LineageConfig
import nextflow.lineage.fs.LinPathFactory
import nextflow.lineage.model.v1beta1.Checksum
import nextflow.lineage.model.v1beta1.DataPath
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.Parameter
import nextflow.lineage.model.v1beta1.Workflow
import nextflow.lineage.model.v1beta1.WorkflowRun
import spock.lang.Specification
import spock.lang.TempDir

import static nextflow.lineage.fs.LinPath.*

/**
 * Lineage channel extensions tests
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class LinExtensionImplTest extends Specification {

    @TempDir
    Path tempDir

    Path storeLocation
    Map configMap

    def setup() {
        storeLocation = tempDir.resolve("store")
        configMap = [linage: [enabled: true, store: [location: storeLocation.toString()]]]
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
        def value2 = new FileOutput("/path/tp/file1", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", "taskid", 1234, time, time, ["value1","value2"])
        def key3 = "testKey3"
        def value3 = new FileOutput("/path/tp/file2", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, ["value2", "value3"])
        def key4 = "testKey4"
        def value4 = new FileOutput("/path/tp/file", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", "taskid", 1234, time, time, ["value4","value3"])
        def lidStore = new DefaultLinStore()
        def session = Mock(Session) {
            getConfig() >> configMap
        }
        lidStore.open(LineageConfig.create(session))
        lidStore.save(key, value1)
        lidStore.save(key2, value2)
        lidStore.save(key3, value3)
        lidStore.save(key4, value4)
        def linExt = Spy(new LinExtensionImpl())
        when:
        def results = CH.create()
        linExt.fromLineage(session, results,  [label: ["value2", "value3"]])
        then:
        linExt.getStore(session) >> lidStore
        and:
        results.val == LinPathFactory.create( asUriString(key3) )
        results.val == Channel.STOP

        when:
        results = CH.create()
        linExt.fromLineage(session, results, [taskRun: "taskid", label: ["value4"]])
        then:
        linExt.getStore(session) >> lidStore
        and:
        results.val == LinPathFactory.create( asUriString(key4) )
        results.val == Channel.STOP

        when:
        results = CH.create()
        linExt.fromLineage(session, results, [workflowRun: "testkey", taskRun: "taskid", label: "value2"])
        then:
        linExt.getStore(session) >> lidStore
        and:
        results.val == LinPathFactory.create( asUriString(key2) )
        results.val == Channel.STOP


    }
}
