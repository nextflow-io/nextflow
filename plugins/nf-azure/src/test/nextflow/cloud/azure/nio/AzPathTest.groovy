/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cloud.azure.nio


import java.time.OffsetDateTime

import com.azure.storage.blob.BlobServiceClient
import com.azure.storage.blob.models.BlobItem
import com.azure.storage.blob.models.BlobItemProperties
import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzPathTest extends Specification {

    @Shared
    Map<String,AzFileSystem> cache = new HashMap<>()

    private AzPath azpath(String path) {
        def container = path.tokenize('/')[0]
        def fs = cache.get(container)
        if( !fs ) {
            def provider = Spy(AzFileSystemProvider)
            provider.@env = [
                    (AzFileSystemProvider.AZURE_STORAGE_ACCOUNT_NAME):'foo',
                    (AzFileSystemProvider.AZURE_STORAGE_ACCOUNT_KEY):'12345' ]
            fs = new AzFileSystem(provider, GroovyMock(BlobServiceClient), container)
            cache.put(container, fs)
        }
        return new AzPath(fs, path)
    }

    @Unroll
    def 'should create a path: #objectName'() {

        when:
        def path = azpath(objectName)
        then:
        path.toString() == expected
        path.@directory == dir

        where:
        objectName              | expected              | dir
        '/bucket/file.txt'      | '/bucket/file.txt'    | false
        '/bucket/a/b/c'         | '/bucket/a/b/c'       | false
        '/bucket/a/b/c/'        | '/bucket/a/b/c'       | true
        '/bucket'               | '/bucket'             | true
        '/bucket/'              | '/bucket'             | true

    }

    @Unroll
    def 'should validate blob constructor and cached attributes'() {
        when:
        def path = azpath(PATH)
        // pre-populate the attributes cache so the slashless case (no directory marker
        // on the path itself) doesn't have to resolve a real blob over the network
        path.setAttributes(new AzFileAttributes(directory: IS_DIRECTORY))
        then:
        path.containerName == CONTAINER
//        path.objectName == BLOB
        path.isDirectory() == IS_DIRECTORY
        path.isContainer() == IS_CONTAINER
        and:
        path.getContainerName() == path.checkContainerName()
//        path.getObjectName() == path.blobName()

        where:
        PATH                    | CONTAINER     | BLOB          | IS_DIRECTORY  | IS_CONTAINER
        '/alpha/beta/delta'     | 'alpha'       | 'beta/delta'  | false         | false
        '/alpha/beta/delta/'    | 'alpha'       | 'beta/delta/' | true          | false
        '/alpha/'               | 'alpha'       | null          | true          | true
        '/alpha'                | 'alpha'       | null          | true          | true

    }

    def 'should validate equals and hashCode'() {

        when:
        def path1 = azpath('/bucket/some/file-name.txt')
        def path2 = azpath('/bucket/some/file-name.txt')
        def path3 = azpath('/bucket/other/file-name.txt')
        def path4 = azpath('/bucket2/some/file-name.txt')

        then:
        path1 == path2
        path1 != path3
        path1 != path4

        path1.hashCode() == path2.hashCode()
        path1.hashCode() != path3.hashCode()

        when:
        def rel1 = azpath('file.txt')
        def rel2 = azpath('file.txt')
        then:
        rel1 == rel2
        rel1.hashCode() == rel2.hashCode()
    }

    def 'should validate isAbsolute'() {
        when:
        def path1 = azpath('/some/file-name.txt')
        def path2 = azpath('file-name.txt')

        then:
        path1.isAbsolute()
        !path2.isAbsolute()
    }

    def 'should validate getRoot'() {
        when:
        def path1 = azpath('/bucket/some/file-name.txt')
        def path2 = azpath('file-name.txt')

        then:
        path1.root == azpath('/bucket')
        path1.root.toString() == '/bucket'
        path2.root == null
    }

    @Unroll
    def 'should return bucket name and id' () {
        when:
        def p = azpath(path)
        then:
        p.isContainer() == expected
        p.getContainerName() == bucketName

        where:
        path                                    | expected  |  bucketName
        '/nxf-bucket/file-name.txt'             | false     | 'nxf-bucket'
        '/nxf-bucket/some/data/file-name.txt'   | false     | 'nxf-bucket'
        'file-name.txt'                         | false     | null
        '/nxf-bucket'                           | true      | 'nxf-bucket'
    }

    @Unroll
    def 'should validate getFileName'() {
        expect:
        azpath(path).getFileName() == azpath(fileName)

        where:
        path                                    | fileName
        '/nxf-bucket/file-name.txt'             | 'file-name.txt'
        '/nxf-bucket/some/data/file-name.txt'   | 'file-name.txt'
        'file-name.txt'                         | 'file-name.txt'
        '/nxf-bucket'                           | 'nxf-bucket'
    }


    @Unroll
    def 'should validate getParent: #path'() {
        expect:
        azpath(path).getParent() == (parent ? azpath(parent) : null)

        where:
        path                                    | parent
        '/nxf-bucket/some/data/file-name.txt'   | '/nxf-bucket/some/data'
        '/nxf-bucket/file-name.txt'             | '/nxf-bucket'
        'file-name.txt'                         | null
        '/nxf-bucket'                           | null
    }

    @Unroll
    def 'should validate toUri: #uri'() {
        expect:
        azpath(path).toUri() == new URI(uri)
        azpath(path).toUri().scheme == new URI(uri).scheme
        azpath(path).toUri().authority == new URI(uri).authority
        azpath(path).toUri().path == new URI(uri).path

        where:
        path                            | uri
        '/alpha/some/file.txt'          | 'az://alpha/some/file.txt'
        '/alpha/'                       | 'az://alpha'
        '/alpha'                        | 'az://alpha'
        'some-file.txt'                 | 'az:some-file.txt'
    }


    @Unroll
    def 'should validate toString: #path'() {
        expect:
        azpath(path).toString() == str

        where:
        path                    | str
        '/alpha/some/file.txt'  | '/alpha/some/file.txt'
        '/alpha'                | '/alpha'
        '/alpha/'               | '/alpha'
        'some-file.txt'         | 'some-file.txt'
    }


    @Unroll
    def 'should validate resolve: base:=#base; path=#path'() {

        expect:
        azpath(base).resolve(path) == azpath(expected)
        azpath(base).resolve( azpath(path) ) == azpath(expected)

        where:
        base                        | path                          | expected
        '/nxf-bucket/some/path'     | 'file-name.txt'               | '/nxf-bucket/some/path/file-name.txt'
        '/nxf-bucket/data'          | 'path/file-name.txt'          | '/nxf-bucket/data/path/file-name.txt'
        '/bucket/data'              | '/other/file-name.txt'        | '/other/file-name.txt'
        '/nxf-bucket'               | 'some/file-name.txt'          | '/nxf-bucket/some/file-name.txt'
    }

    @Unroll
    def 'should validate subpath: #expected'() {
        expect:
        azpath(path).subpath(from, to) == azpath(expected)
        where:
        path                                | from  | to    | expected
        '/bucket/some/big/data/file.txt'    | 0     | 2     | 'bucket/some'
        '/bucket/some/big/data/file.txt'    | 1     | 2     | 'some'
        '/bucket/some/big/data/file.txt'    | 4     | 5     | 'file.txt'
    }

    @Unroll
    def 'should validate startsWith'() {
        expect:
        azpath(path).startsWith(prefix) == expected
        azpath(path).startsWith(azpath(prefix)) == expected

        where:
        path                            | prefix            | expected
        '/bucket/some/data/file.txt'    | '/bucket/some'    | true
        '/bucket/some/data/file.txt'    | '/bucket/'        | true
        '/bucket/some/data/file.txt'    | '/bucket'         | true
        '/bucket/some/data/file.txt'    | 'file.txt'        | false
        'data/file.txt'                 | 'data'            | true
        'data/file.txt'                 | 'file.txt'        | false
    }

    def 'should validate endsWith'() {
        expect:
        azpath(path).endsWith(suffix) == expected
        azpath(path).endsWith(azpath(suffix)) == expected

        where:
        path                            | suffix            | expected
        '/bucket/some/data/file.txt'    | 'file.txt'        | true
        '/bucket/some/data/file.txt'    | 'data/file.txt'   | true
        '/bucket/some/data/file.txt'    | '/data/file.txt'  | false
        '/bucket/some/data/file.txt'    | '/bucket'         | false
        'data/file.txt'                 | 'data'            | false
        'data/file.txt'                 | 'file.txt'        | true
    }


    def 'should validate normalise'() {
        expect:
        azpath(path).normalize() == azpath(expected)
        where:
        path                            | expected
        '/bucket/some/data/file.txt'    | '/bucket/some/data/file.txt'
        '/bucket/some/../file.txt'      | '/bucket/file.txt'
        'bucket/some/../file.txt'       | 'bucket/file.txt'
        'file.txt'                      | 'file.txt'

    }

    @Unroll
    def 'should validate resolveSibling' () {
        expect:
        azpath(base).resolveSibling(path) == azpath(expected)
        azpath(base).resolveSibling(azpath(path)) == azpath(expected)

        where:
        base                    | path                          | expected
        '/bucket/some/path'     | 'file-name.txt'               | '/bucket/some/file-name.txt'
        '/bucket/data'          | 'path/file-name.txt'          | '/bucket/path/file-name.txt'
        '/bucket/data'          | '/other/file-name.txt'        | '/other/file-name.txt'
        '/bucket'               | 'some/file-name.txt'          | '/some/file-name.txt'
    }

    @Unroll
    def 'should validate relativize' () {
        expect:
        azpath(path).relativize(azpath(other)) == azpath(expected)
        where:
        path                    | other                                 | expected
        '/nxf-bucket/some/path' | '/nxf-bucket/some/path/data/file.txt' | 'data/file.txt'
    }

    def 'should validate toAbsolutePath' () {
        expect:
        azpath('/bucket/data/file.txt').toAbsolutePath() == azpath('/bucket/data/file.txt')

        when:
        azpath('file.txt').toAbsolutePath()
        then:
        thrown(UnsupportedOperationException)
    }

    def 'should validate toRealPath' () {
        expect:
        azpath('/bucket/data/file.txt').toRealPath() == azpath('/bucket/data/file.txt')

        when:
        azpath('file.txt').toRealPath()
        then:
        thrown(UnsupportedOperationException)
    }

    def 'should validate iterator' () {
        given:
        def itr = azpath('/nxf-bucket/some/file-name.txt').iterator()
        expect:
        itr.hasNext()
        itr.next() == azpath('nxf-bucket')
        itr.hasNext()
        itr.next() == azpath('some')
        itr.hasNext()
        itr.next() == azpath('file-name.txt')
        !itr.hasNext()

    }

    def 'isDirectory should reuse cached attributes for a path without a trailing slash'() {
        given: 'a slashless path whose resolved attributes report a directory'
        def path = azpath('/pipeline/output')
        def attrs = new AzFileAttributes(size: 0, objectId: '/pipeline/output', directory: true)
        path.setAttributes(attrs)

        expect:
        !path.@directory        // the trailing-slash field alone says "not a directory"
        path.isDirectory()      // ...but the resolved attributes identify it as a directory
        path.attributesCache().is(attrs) // retain the one-shot cache for a later readAttributes call
    }

    def 'isDirectory should not cache attributes resolved from storage'() {
        given:
        def fs = Mock(AzFileSystem)
        fs.getContainerName() >> 'pipeline'
        def path = new AzPath(fs, '/pipeline/output')

        when:
        def first = path.isDirectory()
        def second = path.isDirectory()

        then:
        first
        second
        2 * fs.readAttributes(path) >> new AzFileAttributes(directory: true)
        0 * _
    }

    def 'isDirectory should return false for a relative path'() {
        given:
        def path = azpath('some/file.txt')

        expect:
        !path.isDirectory()
    }

    def 'should recognise hdi_isfolder markers regardless of reported size'() {
        expect:
        AzFileAttributes.isDirectoryMarker(metadata) == directory

        where:
        metadata                    || directory
        [hdi_isfolder: 'true']      || true
        [hdi_isfolder: 'false']     || false
        [:]                         || false
    }

    def 'should classify a non-empty listed hdi_isfolder marker as a directory'() {
        given:
        def properties = new BlobItemProperties()
                .setContentLength(42L)
                .setCreationTime(OffsetDateTime.now())
                .setLastModified(OffsetDateTime.now())
        def item = new BlobItem()
                .setName('output')
                .setMetadata([hdi_isfolder: 'true'])
                .setProperties(properties)

        when:
        def attrs = new AzFileAttributes('pipeline', item)

        then:
        attrs.isDirectory()
        !attrs.isRegularFile()
        attrs.size() == 0L
    }

}
