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

import java.nio.file.Files
import java.nio.file.Path

import nextflow.data.cid.model.Checksum
import nextflow.data.cid.model.TaskOutput
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
        historyLog instanceof CidHistoryFile
    }

    def "save should store value in the correct file location"() {
        given:
        def key = "testKey"
        def value = new TaskOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "cid://source", 1234)
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
        def value = new TaskOutput("/path/to/file", new Checksum("hash_value", "hash_algorithm", "standard"), "cid://source", 1234)
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
}
