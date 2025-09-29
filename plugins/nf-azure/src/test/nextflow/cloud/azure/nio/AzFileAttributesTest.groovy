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
        given:
        def mockClient = GroovyMock(BlobContainerClient) {
            getBlobContainerName() >> 'test-container'
        }

        when:
        def attrs = new AzFileAttributes(mockClient, blobName)

        then:
        attrs.isDirectory() == expectedDirectory
        attrs.isRegularFile() != expectedDirectory
        attrs.fileKey().endsWith("/$blobName")

        where:
        blobName                    | expectedDirectory | comment
        'normal-file.txt'          | false             | 'Regular file without slash'
        'normal-file'              | false             | 'Regular file without slash'
        'problematic-file.txt/'    | true              | 'Path with trailing slash is directory'
        'directory/'               | true              | 'Path with trailing slash is directory'
        'file.log/'                | true              | 'Path with trailing slash is directory'
        'path/to/file.dat/'        | true              | 'Path with trailing slash is directory'
        '/'                        | true              | 'Root slash is directory'
        'multiple///'              | true              | 'Path ending with slashes is directory'
        'has.extension.txt/'       | true              | 'Path with slash is directory regardless of extension'
        'log.2024-01-01.txt/'      | true              | 'Path with slash is directory regardless of extension'
    }

    def 'should validate directory detection for paths without slash'() {
        given:
        def mockClient = GroovyMock(BlobContainerClient) {
            getBlobContainerName() >> 'my-container'
        }

        when:
        def attrs = new AzFileAttributes(mockClient, 'some-directory-without-slash')

        then:
        attrs.isDirectory() == false
        attrs.isRegularFile()
        attrs.fileKey().endsWith('/some-directory-without-slash')
    }

    def 'should handle edge cases in directory detection'() {
        given:
        def mockClient = GroovyMock(BlobContainerClient) {
            getBlobContainerName() >> 'test-container'
        }

        expect:
        new AzFileAttributes(mockClient, 'regular-file/').isDirectory() == true
        new AzFileAttributes(mockClient, 'file.txt/').isDirectory() == true
        new AzFileAttributes(mockClient, '/').isDirectory() == true
        new AzFileAttributes(mockClient, 'multiple///').isDirectory() == true
        new AzFileAttributes(mockClient, 'no-slash').isDirectory() == false
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