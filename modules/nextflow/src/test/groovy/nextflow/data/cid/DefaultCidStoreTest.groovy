package nextflow.data.cid

import nextflow.data.config.DataConfig
import spock.lang.Specification
import spock.lang.TempDir

import java.nio.file.Files
import java.nio.file.Path
import java.util.function.Consumer

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
        cidStore.getPath() == storeLocation
        cidStore.getMetadataPath() == metaLocation
        historyLog != null
        historyLog instanceof CidHistoryFile
    }

    def "save should store value in the correct file location"() {
        given:
        def key = "testKey"
        def value = "testValue"
        def cidStore = new DefaultCidStore()
        cidStore.open(config)

        when:
        cidStore.save(key, value)

        then:
        def filePath = metaLocation.resolve("$key/.data.json")
        Files.exists(filePath)
        filePath.text == value
    }

    def "load should retrieve stored value correctly"() {
        given:
        def key = "testKey"
        def value = "testValue"
        def cidStore = new DefaultCidStore()
        cidStore.open(config)
        cidStore.save(key, value)

        expect:
        cidStore.load(key) == value
    }

    def "load should return null if key does not exist"() {
        given:
        def cidStore = new DefaultCidStore()
        cidStore.open(config)

        expect:
        cidStore.load("nonexistentKey") == null
    }

    def "list should invoke consumer for each stored key"() {
        given:
        def cidStore = new DefaultCidStore()
        cidStore.open(config)
        cidStore.save("key1", "value1")
        cidStore.save("key2", "value2")

        def collectedKeys = []
        Consumer<String> consumer = { collectedKeys << it }

        when:
        cidStore.list("key1", consumer)

        then:
        collectedKeys.contains("key1/.data.json")
    }
}