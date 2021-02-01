package nextflow.cloud.azure.nio

import java.nio.file.FileSystemAlreadyExistsException

import com.azure.storage.blob.BlobServiceClient
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzFileSystemProviderTest extends Specification {

    def 'should return azure storage scheme'() {
        given:
        def provider = new AzFileSystemProvider()
        expect:
        provider.getScheme() == 'az'
    }

    @Unroll
    def 'should return a azure blob path' () {
        given:
        def fs = Mock(AzFileSystem)
        fs.getContainerName() >> 'bucket'
        def provider = Spy(AzFileSystemProvider)

        when:
        def path = provider.getPath(new URI(uri))
        then:
        1 * provider.getFileSystem0(_, true) >> fs
        path == new AzPath(fs, expected)

        where:
        uri                             | expected
        'az:///'                        | '/'
        'az://bucket'                   | '/bucket/'
        'az://bucket/'                  | '/bucket/'
        'az://bucket/this/and/that'     | '/bucket/this/and/that'
        'az://bucket/this/and/that/'    | '/bucket/this/and/that/'

    }

    @Unroll
    def 'should get a azure storage path' () {
        given:
        def fs = Mock(AzFileSystem)
        fs.getContainerName() >> bucket
        def provider = Spy(AzFileSystemProvider)

        when:
        def path = provider.getPath(objectName)
        then:
        1 * provider.getFileSystem0(bucket, true) >> fs
        path == new AzPath(fs, expected)

        where:
        bucket              | objectName            | expected
        'bucket'            | 'bucket'              | '/bucket/'
        'bucket'            | 'bucket/'             | '/bucket/'
        'bucket'            | 'bucket/a/b'          | '/bucket/a/b'
        'bucket'            | 'bucket/a/b/'         | '/bucket/a/b/'
    }

    def 'should return the bucket given a URI'() {
        given:
        def provider = new AzFileSystemProvider()

        expect:
        provider.getContainerName(new URI('az://bucket/alpha/bravo')) == 'bucket'
        provider.getContainerName(new URI('az://BUCKET/alpha/bravo')) == 'bucket'

        when:
        provider.getContainerName(new URI('s3://xxx'))
        then:
        thrown(IllegalArgumentException)

        when:
        provider.getContainerName(new URI('az:/alpha/bravo'))
        then:
        thrown(IllegalArgumentException)

        when:
        provider.getContainerName(new URI('/alpha/bravo'))
        then:
        thrown(IllegalArgumentException)
    }

    def 'should create a new file system with account key'() {

        given:
        def storage = GroovyMock(BlobServiceClient)
        def provider = Spy(AzFileSystemProvider)
        and:
        def NAME = 'xyz'
        def KEY = '1234'
        and:
        def uri = new URI('az://bucket-example/alpha/bravo')
        def env = [AZURE_STORAGE_ACCOUNT_NAME: NAME, AZURE_STORAGE_ACCOUNT_KEY: KEY]

        when:
        def fs = provider.newFileSystem(uri, env)
        then:
        1 * provider.createBlobServiceWithKey(NAME, KEY) >> storage
        fs.containerName == 'bucket-example'
        fs.provider() == provider
        provider.getFileSystem(uri) == fs

        when:
        provider.newFileSystem(uri, env)
        then:
        thrown(FileSystemAlreadyExistsException)
    }

    def 'should create a new file system with sas token'() {

        given:
        def storage = GroovyMock(BlobServiceClient)
        def provider = Spy(AzFileSystemProvider)
        and:
        def NAME = 'xyz'
        def TOKEN = '1234'
        and:
        def uri = new URI('az://bucket-example/alpha/bravo')
        def env = [AZURE_STORAGE_ACCOUNT_NAME: NAME, AZURE_STORAGE_SAS_TOKEN: TOKEN]

        when:
        def fs = provider.newFileSystem(uri, env)
        then:
        1 * provider.createBlobServiceWithToken(NAME, TOKEN) >> storage
        fs.containerName == 'bucket-example'
        fs.provider() == provider
        provider.getFileSystem(uri) == fs

        when:
        provider.newFileSystem(uri, env)
        then:
        thrown(FileSystemAlreadyExistsException)
    }

    def 'should create directory' () {

        given:
        def storage = GroovyMock(BlobServiceClient)
        and:
        def NAME = 'xyz'
        def KEY = '1234'
        def CONTAINER = 'foo'
        and:
        def config = [AZURE_STORAGE_ACCOUNT_NAME: NAME, AZURE_STORAGE_ACCOUNT_KEY: KEY]
        def provider = Spy(AzFileSystemProvider)
        def FS = Mock(AzFileSystem)

        when:
        def fs= provider.newFileSystem0(CONTAINER, config)
        then:
        1 * provider.createBlobServiceWithKey(NAME, KEY) >> storage
        1 * provider.createFileSystem(_, CONTAINER, config) >> FS
        and:
        fs == FS

    }
    
}
