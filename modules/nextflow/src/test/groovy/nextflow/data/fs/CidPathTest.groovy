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

package nextflow.data.fs

import java.nio.file.Path

import spock.lang.Shared
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CidPathTest extends Specification {

    @Shared def BASE = Path.of('/some/base/data')
    @Shared def fs = Mock(CidFileSystem) { getBasePath()>>BASE }

    def 'should create cid path' () {
        when:
        def cid1 = new CidPath(fs, '1234567890')
        then:
        cid1.getFileSystem() == fs
        cid1.toRealPath() == BASE.resolve('1234567890')
        cid1.getTargetSystem() == BASE.fileSystem
        cid1.fileId == new CidFileId('1234567890')

        when:
        def cid2 = new CidPath(fs, '1234567890','foo','bar')
        then:
        cid2.getFileSystem() == fs
        cid2.toRealPath() == BASE.resolve('1234567890/foo/bar')
        cid2.getTargetSystem() == BASE.fileSystem
        cid2.fileId == new CidFileId('1234567890/foo/bar')

        when:
        def cid3 = new CidPath(fs, '/' )
        then:
        cid3.getFileSystem() == fs
        cid3.toRealPath() == BASE
        cid3.getTargetSystem() == BASE.fileSystem
        cid3.fileId == new CidFileId('/')
    }

    def 'should get file name' () {
        when:
        def cid1 = new CidPath(fs, '1234567890/this/file.bam')
        then:
        cid1.getFileName() == new CidPath(null, 'file.bam')
    }

    def 'should get file parent' () {
        when:
        def cid1 = new CidPath(fs, '1234567890/this/file.bam')
        then:
        cid1.getParent() == new CidPath(fs, '1234567890/this')
        cid1.getParent().getParent() == new CidPath(fs, '1234567890')
        cid1.getParent().getParent().getParent() == new CidPath(fs, "/")
        cid1.getParent().getParent().getParent().getParent() == null
    }

    @Unroll
    def 'should get name count' () {
        expect:
        new CidPath(fs, PATH).getNameCount() == EXPECTED
        where:
        PATH        | EXPECTED
        '/'         | 0
        '123'       | 1
        '123/a'     | 2
        '123/a/'    | 2
        '123/a/b'   | 3
    }

    @Unroll
    def 'should get name by index' () {
        expect:
        new CidPath(fs, PATH).getName(INDEX) == EXPECTED
        where:
        PATH        | INDEX | EXPECTED
        '123'       | 0     | new CidPath(fs, '123')
        '123/a'     | 1     | new CidPath(null, 'a')
        '123/a/'    | 1     | new CidPath(null, 'a')
        '123/a/b'   | 2     | new CidPath(null, 'b')
    }

    @Unroll
    def 'should get subpath' () {
        expect:
        new CidPath(fs, PATH).subpath(BEGIN,END) == EXPECTED
        where:
        PATH        | BEGIN | END   | EXPECTED
        '123'       | 0     | 1     | new CidPath(fs, '123')
        '123/a'     | 0     | 2     | new CidPath(fs, '123/a')
        '123/a/'    | 0     | 2     | new CidPath(fs, '123/a')
        '123/a'     | 1     | 2     | new CidPath(null, 'a')
        '123/a/'    | 1     | 2     | new CidPath(null, 'a')
        '123/a/b'   | 2     | 3     | new CidPath(null, 'b')
        '123/a/b'   | 1     | 3     | new CidPath(null, 'a/b')
    }

    def 'should normalize a path' () {
        expect:
        new CidPath(fs, '123').normalize() == new CidPath(fs, '123')
        new CidPath(fs, '123/a/b').normalize() == new CidPath(fs, '123/a/b')
        new CidPath(fs, '123/./a/b').normalize() == new CidPath(fs, '123/a/b')
        new CidPath(fs, '123/a/../a/b').normalize() == new CidPath(fs, '123/a/b')
    }

    @Unroll
    def 'should validate startWith' () {
        expect:
        new CidPath(fs,PATH).startsWith(OTHER) == EXPECTED
        where:
        PATH            | OTHER         | EXPECTED
        '12345/a/b'     | '12345'       | true
        '12345/a/b'     | '12345/a'     | true
        '12345/a/b'     | '12345/a/b'   | true
        and:
        '12345/a/b'     | '12345/b'     | false
        '12345/a/b'     | 'xyz'         | false
    }

    @Unroll
    def 'should validate endsWith' () {
        expect:
        new CidPath(fs,PATH).endsWith(OTHER) == EXPECTED
        where:
        PATH            | OTHER         | EXPECTED
        '12345/a/b'     | 'b'           | true
        '12345/a/b'     | 'a/b'         | true
        '12345/a/b'     | '12345/a/b'   | true
        and:
        '12345/a/b'     | '12345/b'     | false
        '12345/a/b'     | 'xyz'         | false
    }

    def 'should validate isAbsolute' () {
        expect:
        new CidPath(fs,'1234/a/b/c').isAbsolute()
        new CidPath(fs,'1234/a/b/c').getRoot().isAbsolute()
        new CidPath(fs,'1234/a/b/c').getParent().isAbsolute()
        new CidPath(fs,'1234/a/b/c').normalize().isAbsolute()
        new CidPath(fs,'1234/a/b/c').getName(0).isAbsolute()
        new CidPath(fs,'1234/a/b/c').subpath(0,2).isAbsolute()
        and:
        !new CidPath(fs,'1234/a/b/c').getFileName().isAbsolute()
        !new CidPath(fs,'1234/a/b/c').getName(1).isAbsolute()
        !new CidPath(fs,'1234/a/b/c').subpath(1,3).isAbsolute()
    }

    @Unroll
    def 'should get root path' () {
        expect:
        new CidPath(fs,PATH).getRoot() == new CidPath(fs,EXPECTED)
        where:
        PATH            | EXPECTED
        '12345'         | '/'
        '12345/a'       | '/'
    }

    def 'should resolve path' () {
        when:
        def cid1 = new CidPath(fs, '123/a/b/c')
        def cid2 = new CidPath(fs, '321/x/y/z')
        def rel1 = new CidPath(null, 'foo')
        def rel2 = new CidPath(null, 'bar/')
        
        then:
        cid1.resolve(cid2) == cid2
        cid2.resolve(cid1) == cid1
        and:
        cid1.resolve(rel1) == new CidPath(fs,'123/a/b/c/foo')
        cid1.resolve(rel2) == new CidPath(fs,'123/a/b/c/bar')
        and:
        rel1.resolve(rel2) == new CidPath(null, 'foo/bar')
        rel2.resolve(rel1) == new CidPath(null, 'bar/foo')
    }

    def 'should resolve path as string' () {
        given:
        def pr = Mock(CidFileSystemProvider)
        def fs = Mock(CidFileSystem) { getBasePath()>>BASE; provider()>>pr }
        def cid1 = new CidPath(fs, '123/a/b/c')

        expect:
        cid1.resolve('x/y') == new CidPath(fs, '123/a/b/c/x/y')
        cid1.resolve('/x/y/') == new CidPath(fs, '123/a/b/c/x/y')

        when:
        def result = cid1.resolve('cid://321')
        then:
        pr.getPath(CidPath.asUri('cid://321')) >> new CidPath(fs, '321')
        and:
        result == new CidPath(fs, '321')
    }
    
    @Unroll
    def 'should get to uri string' () {
        expect:
        new CidPath(fs, PATH).toUriString() == EXPECTED
        where:
        PATH            | EXPECTED
        '/'             | 'cid:///'
        '/1234'         | 'cid://1234'
        '/1234/a/b/c'   | 'cid://1234/a/b/c'
    }
}
