/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.data.cid

import nextflow.data.cid.model.DataPath
import nextflow.data.cid.model.DataOutput
import nextflow.data.cid.model.Parameter
import nextflow.data.cid.model.Workflow
import nextflow.data.cid.model.WorkflowRun

import java.nio.file.Files
import java.nio.file.Path
import java.time.Instant

import nextflow.data.cid.model.Checksum
import nextflow.data.cid.serde.CidEncoder
import nextflow.data.config.DataConfig
import spock.lang.Specification
import spock.lang.TempDir
/**
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class DefaultCidStoreTest extends Specification {

    @TempDir
    Path tempDir

    Path storeLocation
    Path metaLocation
    DataConfig config

    def setup() {
        storeLocation = tempDir.resolve("store")
        metaLocation = storeLocation.resolve(".meta")
        def configMap = [enabled: true, store: [location: storeLocation.toString(), logLocation: storeLocation.resolve(".log").toString()]]
        config = new DataConfig(configMap)
    }

    def 'should open store'() {
        given:
        def cidStore = new DefaultCidStore()
        when:
        cidStore.open(config)
        def historyLog = cidStore.getHistoryLog()
        then:
        cidStore.getMetadataPath() == metaLocation
        historyLog != null
        historyLog instanceof DefaultCidHistoryLog
    }

    def "save should store value in the correct file location"() {
        given:
        def key = "testKey"
        def value = new DataOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "cid://source", "cid://workflow", "cid://task", 1234)
        def cidStore = new DefaultCidStore()
        cidStore.open(config)

        when:
        cidStore.save(key, value)

        then:
        def filePath = metaLocation.resolve("$key/.data.json")
        Files.exists(filePath)
        filePath.text == new CidEncoder().encode(value)
    }

    def "load should retrieve stored value correctly"() {
        given:
        def key = "testKey"
        def value = new DataOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "cid://source", "cid://workflow", "cid://task", 1234)
        def cidStore = new DefaultCidStore()
        cidStore.open(config)
        cidStore.save(key, value)

        expect:
        cidStore.load(key).toString() == value.toString()
    }

    def "load should return null if key does not exist"() {
        given:
        def cidStore = new DefaultCidStore()
        cidStore.open(config)

        expect:
        cidStore.load("nonexistentKey") == null
    }

    def 'should query' () {
        given:
        def uniqueId = UUID.randomUUID()
        def time = Instant.ofEpochMilli(1234567)
        def mainScript = new DataPath("file://path/to/main.nf", new Checksum("78910", "nextflow", "standard"))
        def workflow = new Workflow([mainScript],"https://nextflow.io/nf-test/", "123456" )
        def key = "testKey"
        def value1 = new WorkflowRun(workflow, uniqueId.toString(), "test_run", [ new Parameter("String", "param1", "value1"), new Parameter("String", "param2", "value2")] )
        def key2 = "testKey2"
        def value2 = new DataOutput("/path/tp/file1", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, [key1:"value1", key2:"value2"])
        def key3 = "testKey3"
        def value3 = new DataOutput("/path/tp/file2", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, [key2:"value2", key3:"value3"])
        def key4 = "testKey4"
        def value4 = new DataOutput("/path/tp/file", new Checksum("78910", "nextflow", "standard"), "testkey", "testkey", null, 1234, time, time, [key3:"value3", key4:"value4"])

        def cidStore = new DefaultCidStore()
        cidStore.open(config)
        cidStore.save(key, value1)
        cidStore.save(key2, value2)
        cidStore.save(key3, value3)
        cidStore.save(key4, value4)

        when:
        def results3 = cidStore.search("type=DataOutput&annotations.key2=value2")
        then:
        results3.size() == 2
    }


}
