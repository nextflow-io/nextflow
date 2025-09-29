package nextflow.cloud.azure.nio

import com.azure.storage.blob.BlobContainerClient
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Unit tests for AzFileAttributes class
 */
class AzFileAttributesTest extends Specification {

    def 'should create root attributes correctly'() {
        when:
        def attrs = AzFileAttributes.root()

        then:
        attrs.isDirectory()
        !attrs.isRegularFile()
        attrs.size() == 0
        attrs.fileKey() == '/'
    }

    @Unroll
    def 'should validate directory detection with blobName: #blobName'() {
        when:
        def attrs = new AzFileAttributes()
        attrs.objectId = "/test-container/$blobName"
        attrs.directory = expectedDirectory
        attrs.size = expectedDirectory ? 0 : 100

        then:
        attrs.isDirectory() == expectedDirectory
        attrs.isRegularFile() != expectedDirectory
        attrs.fileKey().endsWith("/$blobName")

        where:
        blobName                    | expectedDirectory | comment
        'problematic-file.txt/'    | true              | 'Path with trailing slash is directory'
        'directory/'               | true              | 'Path with trailing slash is directory'
        'file.log/'                | true              | 'Path with trailing slash is directory'
        'path/to/file.dat/'        | true              | 'Path with trailing slash is directory'
        '/'                        | true              | 'Root slash is directory'
        'multiple///'              | true              | 'Path ending with slashes is directory'
        'has.extension.txt/'       | true              | 'Path with slash is directory regardless of extension'
        'log.2024-01-01.txt/'      | true              | 'Path with slash is directory regardless of extension'
    }



    def 'should verify equality and hashCode methods work correctly'() {
        given:
        def attrs1 = AzFileAttributes.root()
        def attrs2 = AzFileAttributes.root()

        when:
        def equals1 = attrs1.equals(attrs2)
        def hash1 = attrs1.hashCode()
        def hash2 = attrs2.hashCode()

        then:
        equals1 == true
        hash1 == hash2
    }
}