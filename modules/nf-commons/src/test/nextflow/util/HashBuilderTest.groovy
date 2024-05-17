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
 */

package nextflow.util

import java.nio.file.Files
import java.nio.file.Paths

import com.google.common.hash.Hashing
import nextflow.Global
import nextflow.Session
import org.apache.commons.codec.digest.DigestUtils
import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class HashBuilderTest extends Specification {


    def testHashContent() {
        setup:
        def path1 = Files.createTempFile('test-hash-content',null)
        def path2 = Files.createTempFile('test-hash-content',null)
        def path3 = Files.createTempFile('test-hash-content',null)

        path1.text = '''
            line 1
            line 2
            line 3 the file content
            '''


        path2.text = '''
            line 1
            line 2
            line 3 the file content
            '''

        path3.text = '''
            line 1
            line 1
            line 1 the file content
            '''

        expect:
        HashBuilder.hashContent(path1) == HashBuilder.hashContent(path2)
        HashBuilder.hashContent(path1) != HashBuilder.hashContent(path3)
        HashBuilder.hashContent(path1, Hashing.md5()) == HashBuilder.hashContent(path2,Hashing.md5())
        HashBuilder.hashContent(path1, Hashing.md5()) != HashBuilder.hashContent(path3,Hashing.md5())

        cleanup:
        path1.delete()
        path2.delete()
        path3.delete()

    }
    
    def 'should validate is asset file'() {
        when:
        def BASE = Paths.get("/some/pipeline/dir")
        and:
        Global.session = Mock(Session) { getBaseDir() >> BASE }
        then:
        !HashBuilder.isAssetFile(BASE.resolve('foo'))


        when:
        Global.session = Mock(Session) {
            getBaseDir() >> BASE
            getCommitId() >> '123456'
        }
        then:
        HashBuilder.isAssetFile(BASE.resolve('foo'))
        and:
        !HashBuilder.isAssetFile(Paths.get('/other/dir'))
    }

    def 'should hash file content'() {
        given:
        def EXPECTED = '64ec88ca00b268e5ba1a35678a1b5316d212f4f366b2477232534a8aeca37f3c'
        def file = TestHelper.createInMemTempFile('foo', 'Hello world')
        expect:
        HashBuilder.hashFileSha256Impl0(file) == EXPECTED
        and:
        HashBuilder.hashFileSha256Impl0(file) == DigestUtils.sha256Hex(file.bytes)
    }

    def 'should hash dir content with sha256'() {
        given:
        def folder = TestHelper.createInMemTempDir()
        folder.resolve('dir1').mkdir()
        folder.resolve('dir2').mkdir()
        and:
        folder.resolve('dir1/foo').text = "I'm foo"
        folder.resolve('dir1/bar').text = "I'm bar"
        folder.resolve('dir1/xxx/yyy').mkdirs()
        folder.resolve('dir1/xxx/foo1').text = "I'm foo within xxx"
        folder.resolve('dir1/xxx/yyy/bar1').text = "I'm bar within yyy"
        and:
        folder.resolve('dir2/foo').text = "I'm foo"
        folder.resolve('dir2/bar').text = "I'm bar"
        folder.resolve('dir2/xxx/yyy').mkdirs()
        folder.resolve('dir2/xxx/foo1').text = "I'm foo within xxx"
        folder.resolve('dir2/xxx/yyy/bar1').text = "I'm bar within yyy"

        when:
        def hash1 = HashBuilder.hashDirSha256(HashBuilder.defaultHasher(), folder.resolve('dir1'), folder.resolve('dir1'))
        and:
        def hash2 = HashBuilder.hashDirSha256(HashBuilder.defaultHasher(), folder.resolve('dir2'), folder.resolve('dir2'))

        then:
        hash1.hash() == hash2.hash()

    }
}
