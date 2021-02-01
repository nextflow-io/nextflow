package nextflow.cloud.azure.nio

import java.nio.file.Paths

import com.azure.storage.blob.BlobServiceClient
import spock.lang.Ignore
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzFileSystemTest extends Specification {


    def 'should get a path' () {
        given:
        final storage = GroovyMock(BlobServiceClient)
        final provider = Spy(AzFileSystemProvider)
        final foo_fs = new AzFileSystem(provider, storage, 'foo')
        final bar_fs = new AzFileSystem(provider, storage, 'bar')

        when:
        def result = foo_fs.getPath(path, more as String[])
        then:
        call * provider.getFileSystem0(_,true) >> { it, create -> it=='foo' ? foo_fs : bar_fs }
        result.toUriString()

        where:
        call| path                  | more          | expected
        1   | '/foo'                | null          | 'az://foo/'
        1   | '/foo/'               | null          | 'az://foo/'
        1   | '/foo/alpha/bravo'    | null          | 'az://foo/alpha/bravo'
        1   | '/foo/alpha/bravo/'   | null          | 'az://foo/alpha/bravo/'
        1   | '/bar'                | null          | 'az://bar/'
        1   | '/bar'                | ['a','b']     | 'az://bar/a/b'
        1   | '/bar/'               | ['a/','b/']   | 'az://bar/a/b'
        1   | '/bar/'               | ['/a','/b']   | 'az://bar/a/b'
        0   | 'this/and/that'       | null          | 'az:/this/and/that'
        0   | 'this/and/that'       | 'x/y'         | 'az:/this/and/that/x/y'

    }

    def 'should return root path' () {

        given:
        def provider = Mock(AzFileSystemProvider)
        final storage = GroovyMock(BlobServiceClient)
        def fs = new AzFileSystem(provider, storage, 'bucket-example')

        expect:
        fs.getRootDirectories() == [ new AzPath(fs, '/bucket-example/') ]
    }

    @Ignore
    def 'should return root paths' () {

        given:
        def page = Mock(Page)
        final storage = GroovyMock(BlobServiceClient)
        and:
        def bucket1 = Mock(Bucket); bucket1.getName() >> 'alpha'
        def bucket2 = Mock(Bucket); bucket2.getName() >> 'beta'
        def bucket3 = Mock(Bucket); bucket3.getName() >> 'delta'
        and:
        def provider = Mock(AzFileSystemProvider)
        and:
        def fs = new AzFileSystem(provider, storage, '/')

        when:
        def roots = fs.getRootDirectories()
        then:
        1 * storage.list() >> page
        1 * page.iterateAll() >> { [bucket1, bucket2, bucket3].iterator() }
        3 * provider.getPath(_ as String) >>  { String it -> new AzPath(fs:Mock(AzFileSystem), path: Paths.get("/$it")) }
        roots.size() == 3
        roots[0].toString() == '/alpha'
        roots[1].toString() == '/beta'
        roots[2].toString() == '/delta'
    }

    def 'should test basic properties' () {

        given:
        def BUCKET_NAME = 'bucket'
        def provider = Stub(AzFileSystemProvider)
        final storage = GroovyMock(BlobServiceClient)

        when:
        def fs = new AzFileSystem(provider, storage, BUCKET_NAME)
        then:
        fs.getSeparator() == '/'
        fs.isOpen()
        fs.provider() == provider
        fs.containerName == BUCKET_NAME
        !fs.isReadOnly()
        fs.supportedFileAttributeViews() == ['basic'] as Set

        when:
        fs = new AzFileSystem(provider, storage, '/')
        then:
        fs.containerName == '/'
        fs.isReadOnly()
        fs.isOpen()
    }

    def 'should test getPath' () {
        given:
        def BUCKET_NAME = 'bucket'
        def provider = Stub(AzFileSystemProvider)
        final storage = GroovyMock(BlobServiceClient)
        and:
        def fs = new AzFileSystem(provider, storage, BUCKET_NAME)

        expect:
        fs.getPath('file-name.txt') == new AzPath(fs, Paths.get('file-name.txt'), false)
        fs.getPath('alpha/bravo') == new AzPath(fs, Paths.get('alpha/bravo'), false)
        fs.getPath('/alpha/bravo') == new AzPath(fs, Paths.get('/bucket/alpha/bravo'), false)
        fs.getPath('/alpha','/gamma','/delta') == new AzPath(fs, Paths.get('/bucket/alpha/gamma/delta'), false)
        fs.getPath('/alpha','gamma//','delta//') == new AzPath(fs, Paths.get('/bucket/alpha/gamma/delta'), false)
    }
    
}
