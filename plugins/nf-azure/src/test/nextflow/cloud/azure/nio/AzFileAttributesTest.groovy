package nextflow.cloud.azure.nio

import com.azure.storage.blob.BlobContainerClient
import spock.lang.Specification
import spock.lang.Unroll

/**
 * Unit tests for AzFileAttributes class
 *
 * This focuses on testing the problematic constructor that determines directory status
 * based solely on trailing slashes, which is the bug reported in issue #6427.
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

    /**
     * Tests the problematic constructor: AzFileAttributes(BlobContainerClient client, String blobName)
     * This constructor at line 114 in AzFileAttributes.groovy has the bug:
     * directory = blobName.endsWith('/')
     */
    @Unroll
    def 'should demonstrate trailing slash bug in constructor with blobName: #blobName'() {
        given:
        def mockClient = GroovyMock(BlobContainerClient) {
            getBlobContainerName() >> 'test-container'
        }

        when:
        def attrs = new AzFileAttributes(mockClient, blobName)

        then:
        // Document the current buggy behavior - any path ending with '/' is treated as directory
        attrs.isDirectory() == expectedDirectory
        attrs.isRegularFile() != expectedDirectory
        // Note: fileKey shows "/null/$blobName" due to mocking issues, but directory logic still works
        attrs.fileKey().endsWith("/$blobName")

        where:
        blobName                    | expectedDirectory | comment
        'normal-file.txt'          | false             | 'Correct: regular file without slash'
        'normal-file'              | false             | 'Correct: regular file without slash'
        'problematic-file.txt/'    | true              | 'BUG: Should be false - file with trailing slash'
        'directory/'               | true              | 'BUG: Should require metadata to confirm directory'
        'file.log/'                | true              | 'BUG: Should be false - file with trailing slash'
        'path/to/file.dat/'        | true              | 'BUG: Should be false - file with trailing slash'
        '/'                        | true              | 'BUG: Should require metadata to confirm directory'
        'multiple///'              | true              | 'BUG: Should be false - file with multiple slashes'
        'has.extension.txt/'       | true              | 'BUG: Should be false - clearly a file with extension'
        'log.2024-01-01.txt/'      | true              | 'BUG: Should be false - log file with timestamp'
    }

    def 'should demonstrate the specific false positive case from GitHub issue #6427'() {
        given:
        def mockClient = GroovyMock(BlobContainerClient) {
            getBlobContainerName() >> 'my-container'
        }

        when:
        // This represents the exact scenario from the GitHub issue:
        // A blob path that ends with '/' but is actually a file, not a directory
        def attrs = new AzFileAttributes(mockClient, 'some-problematic-filename/')

        then:
        // CURRENT BEHAVIOR (BUG):
        // The constructor incorrectly treats this as a directory solely because it ends with '/'
        attrs.isDirectory() == true
        !attrs.isRegularFile()
        // Note: fileKey shows "/null/..." due to mocking, but the directory detection bug still exists
        attrs.fileKey().endsWith('/some-problematic-filename/')

        // EXPECTED BEHAVIOR:
        // Should be false unless there's proper metadata (hdi_isfolder='true') or isPrefix=true
        // This test documents the bug that needs to be fixed
    }

    def 'should handle edge cases with the problematic constructor'() {
        given:
        def mockClient = GroovyMock(BlobContainerClient) {
            getBlobContainerName() >> 'test-container'
        }

        expect:
        // Test various edge cases that all demonstrate the trailing slash bug
        new AzFileAttributes(mockClient, 'regular-file/').isDirectory() == true     // BUG
        new AzFileAttributes(mockClient, 'file.txt/').isDirectory() == true         // BUG
        new AzFileAttributes(mockClient, '/').isDirectory() == true                 // BUG
        new AzFileAttributes(mockClient, 'multiple///').isDirectory() == true       // BUG
        new AzFileAttributes(mockClient, 'no-slash').isDirectory() == false         // CORRECT

        // All the above cases ending with '/' are incorrectly treated as directories
        // Only the last case without trailing slash correctly returns false
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