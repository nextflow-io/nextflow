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

import java.nio.file.Files
import java.nio.file.Path
import java.time.Instant
import java.time.OffsetDateTime
import java.time.ZoneOffset

import nextflow.lineage.model.v1beta1.Checksum
import nextflow.lineage.model.v1beta1.DataPath
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.Parameter
import nextflow.lineage.model.v1beta1.Workflow
import nextflow.lineage.model.v1beta1.WorkflowRun
import nextflow.lineage.serde.LinEncoder
import nextflow.lineage.config.LineageConfig
import spock.lang.Specification
import spock.lang.TempDir

/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class DefaultLinStoreTest extends Specification {

    @TempDir
    Path tempDir

    Path storeLocation
    LineageConfig config

    def setup() {
        storeLocation = tempDir.resolve("store")
        def configMap = [enabled: true, store: [location: storeLocation.toString(), logLocation: storeLocation.resolve(".log").toString()]]
        config = new LineageConfig(configMap)
    }

    def 'should open store'() {
        given:
        def store = new DefaultLinStore()
        when:
        store.open(config)
        def historyLog = store.getHistoryLog()
        then:
        store.getLocation() == storeLocation
        historyLog != null
        historyLog instanceof DefaultLinHistoryLog
    }

    def "save should store value in the correct file location"() {
        given:
        def key = "testKey"
        def value = new FileOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "lid://source", "lid://workflow", "lid://task", 1234)
        def lidStore = new DefaultLinStore()
        lidStore.open(config)

        when:
        lidStore.save(key, value)

        then:
        def filePath = storeLocation.resolve("$key/.data.json")
        Files.exists(filePath)
        filePath.text == new LinEncoder().encode(value)
    }

    def "load should retrieve stored value correctly"() {
        given:
        def key = "testKey"
        def value = new FileOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "lid://source", "lid://workflow", "lid://task", 1234)
        def lidStore = new DefaultLinStore()
        lidStore.open(config)
        lidStore.save(key, value)

        expect:
        lidStore.load(key).toString() == value.toString()
    }

    def "load should return null if key does not exist"() {
        given:
        def lidStore = new DefaultLinStore()
        lidStore.open(config)

        expect:
        lidStore.load("nonexistentKey") == null
    }

    def 'should query' () {
        given:
        def uniqueId = UUID.randomUUID()
        def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(1234567), ZoneOffset.UTC)
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript],"https://nextflow.io/nf-test/", "123456" )
        def key = "testKey"
        def value1 = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [ new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")] )
        def key2 = "testKey2"
        def value2 = new FileOutput("/path/tp/file1", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, ["value1", "value2"])
        def key3 = "testKey3"
        def value3 = new FileOutput("/path/tp/file2", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, ["value2", "value3"])
        def key4 = "testKey4"
        def value4 = new FileOutput("/path/tp/file", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, ["value4", "value3"])

        def lidStore = new DefaultLinStore()
        lidStore.open(config)
        lidStore.save(key, value1)
        lidStore.save(key2, value2)
        lidStore.save(key3, value3)
        lidStore.save(key4, value4)

        when:
        def results = lidStore.search( [type:['FileOutput'], labels:['value2']]).toList()
        then:
        results.size() == 2
        results.containsAll([key2,key3])
    }

    def 'should search subkeys' () {
         given:
         def uniqueId = UUID.randomUUID()
         def time = OffsetDateTime.ofInstant(Instant.ofEpochMilli(1234567), ZoneOffset.UTC)
         def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
         def workflow = new Workflow([mainScript], "https://nextflow.io/nf-test/", "123456")
         def key = "testKey"
         def value1 = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")])
         def key2 = "testKey/file1"
         def value2 = new FileOutput("/path/file1", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, ["value1", "value2"])
         def key3 = "testKey/subfolder/file3"
         def value3 = new FileOutput("/path//file2", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, ["value2", "value3"])
         def key4 = "testKey2/file2"
         def value4 = new FileOutput("/path/tp/file", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, ["value4", "value3"])

         def lidStore = new DefaultLinStore()
         lidStore.open(config)
         lidStore.save(key, value1)
         lidStore.save(key2, value2)
         lidStore.save(key3, value3)
         lidStore.save(key4, value4)

         when:
         def results = lidStore.getSubKeys("testKey").toList()
         then:
         results.size() == 2
         results.containsAll([key2, key3])
     }

}
