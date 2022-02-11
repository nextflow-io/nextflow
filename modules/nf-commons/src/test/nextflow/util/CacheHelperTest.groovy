/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
import java.nio.file.attribute.FileTime

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
class CacheHelperTest extends Specification {

    def 'test hashCode' () {

        setup:
        def anInteger = Hashing.murmur3_128().newHasher().putInt(1).hash()
        def aLong = Hashing.murmur3_128().newHasher().putLong(100L).hash()
        def aString = Hashing.murmur3_128().newHasher().putUnencodedChars('Hola').hash()
        def aDouble = Hashing.murmur3_128().newHasher().putDouble(10.5).hash()
        def aBool = Hashing.murmur3_128().newHasher().putBoolean(true).hash()
        def aByteArray = Hashing.murmur3_128().newHasher().putBytes([0x1,0x2,0x3] as byte[]).hash()
        def anObjectArray = Hashing.murmur3_128().newHasher().putInt(1).putInt(2).putInt(3).hash()
        def aMap =  Hashing.murmur3_128().newHasher().putInt(1).putUnencodedChars('String1').putBoolean(true).hash()
        def aList = Hashing.murmur3_128().newHasher().putUnencodedChars('A').putUnencodedChars('B').putUnencodedChars('C').hash()

        def file = Files.createTempFile('test', null)
        file.deleteOnExit()
        file.text = 'Hello'

        def aFile = Hashing.murmur3_128()
                        .newHasher()
                        .putUnencodedChars(file.toAbsolutePath().toString())
                        .putLong(Files.size(file))
                        .putLong(Files.getLastModifiedTime(file).toMillis())
                        .hash()

        expect:
        CacheHelper.hasher(1).hash() == anInteger
        CacheHelper.hasher(100L).hash() == aLong
        CacheHelper.hasher('Hola').hash() == aString
        CacheHelper.hasher(10.5D).hash() == aDouble
        CacheHelper.hasher(true).hash() == aBool
        CacheHelper.hasher([0x1,0x2,0x3] as byte[]).hash() == aByteArray
        CacheHelper.hasher([1,2,3] as Object[]).hash() == anObjectArray
        CacheHelper.hasher( [f1: 1, f2: 'String1', f3: true] ) .hash() == aMap
        CacheHelper.hasher( ['A','B','C'] ).hash() == aList
        CacheHelper.hasher(file).hash() == aFile

        CacheHelper.hasher(['abc',123]).hash() == CacheHelper.hasher(['abc',123]).hash()
        CacheHelper.hasher(['abc',123]).hash() != CacheHelper.hasher([123,'abc']).hash()

        // set hash-code is order invariant
        CacheHelper.hasher(['abc',123] as Set ).hash() == CacheHelper.hasher(['abc',123] as Set ).hash()
        CacheHelper.hasher(['abc',123] as Set ).hash() == CacheHelper.hasher([123,'abc'] as Set ).hash()

        CacheHelper.hasher(['abc',123] as Set ).hash() == CacheHelper.hasher([123,'abc'] as Set ).hash()


    }


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
        CacheHelper.hashContent(path1) == CacheHelper.hashContent(path2)
        CacheHelper.hashContent(path1) != CacheHelper.hashContent(path3)
        CacheHelper.hashContent(path1, Hashing.md5()) == CacheHelper.hashContent(path2,Hashing.md5())
        CacheHelper.hashContent(path1, Hashing.md5()) != CacheHelper.hashContent(path3,Hashing.md5())

        cleanup:
        path1.delete()
        path2.delete()
        path3.delete()

    }


    def testHashOrder () {

        when:
        def h1 = Hashing.murmur3_128().newHasher()
        def h2 = Hashing.murmur3_128().newHasher()
        def c1 = h1.putInt(1).putInt(2).hash()
        def c2 = h2.putInt(2).putInt(1).hash()
        then:
        c1 != c2

        when:
        h1 = Hashing.murmur3_128().newHasher()
        h2 = Hashing.murmur3_128().newHasher()

        c1 = h1.putInt(1).putInt(2).hash()
        h2.putInt(1)
        h2.putInt(2)
        c2 = h2.hash()

        then:
        c1 == c2

    }


    def testSetAndMap() {

        expect:
        [x:1, y:2, z:3 ].hashCode() == [z:3, y:2, x:1 ].hashCode()
        new HashSet(['a','b','z']).hashCode() == new HashSet(['z','b','a']).hashCode()

    }

    def 'should hash file in standard mode' () {

        given:
        def MODE = CacheHelper.HashMode.STANDARD
        def folder = Files.createTempDirectory('test')

        def file = folder.resolve('file.txt'); file.text = 'hello world'
        def hash = CacheHelper .hasher(file, MODE) .hash() .toString()

        // modified modified time does not affect hashcode
        when:
        def time = FileTime.fromMillis(System.currentTimeMillis()-100_000)
        Files.setLastModifiedTime(file, time)
        then:
        CacheHelper .hasher(file, MODE) .hash() .toString() != hash

        cleanup:
        folder?.deleteDir()
    }

    def 'should hash file in lenient mode' () {

        given:
        def MODE = CacheHelper.HashMode.LENIENT
        def folder = Files.createTempDirectory('test')

        def file = folder.resolve('file.txt'); file.text = 'hello world'
        def hash = CacheHelper .hasher(file, MODE) .hash() .toString()

        // modified modified time does not affect hashcode
        when:
        def time = FileTime.fromMillis(System.currentTimeMillis()-100_000)
        Files.setLastModifiedTime(file, time)
        then:
        CacheHelper .hasher(file, MODE) .hash() .toString() == hash

        cleanup:
        folder?.deleteDir()
    }


    def 'should hash file in deep mode' () {

        given:
        def MODE = CacheHelper.HashMode.DEEP
        def folder = Files.createTempDirectory('test')

        when:
        def file = folder.resolve('file.txt'); file.text = 'hello world'
        then:
        CacheHelper .hasher(file, MODE) .hash() .toString() == '0e617feb46603f53b163eb607d4697ab'

        // modified modified time does not affect hashcode
        when:
        def time = FileTime.fromMillis(System.currentTimeMillis()-100_000)
        Files.setLastModifiedTime(file, time)
        then:
        CacheHelper .hasher(file, MODE) .hash() .toString() == '0e617feb46603f53b163eb607d4697ab'

        cleanup:
        folder?.deleteDir()
    }

    def "should get hash mode from value=#VALUE" () {

        expect:
        CacheHelper.HashMode.of(VALUE) == EXPECTED

        where:
        VALUE      | EXPECTED
        null       | null
        true       | null
        'true'     | null
        'false'    | null
        'standard' | CacheHelper.HashMode.STANDARD
        'deep'     | CacheHelper.HashMode.DEEP
        'lenient'  | CacheHelper.HashMode.LENIENT
        'sha256'   | CacheHelper.HashMode.SHA256
    }
    
    def 'should validate is asset file'() {
        when:
        def BASE = Paths.get("/some/pipeline/dir")
        and:
        Global.session = Mock(Session) { getBaseDir() >> BASE }
        then:
        !CacheHelper.isAssetFile(BASE.resolve('foo'))


        when:
        Global.session = Mock(Session) {
            getBaseDir() >> BASE
            getCommitId() >> '123456'
        }
        then:
        CacheHelper.isAssetFile(BASE.resolve('foo'))
        and:
        !CacheHelper.isAssetFile(Paths.get('/other/dir'))
    }


    def 'should hash file content'() {
        given:
        def EXPECTED = '64ec88ca00b268e5ba1a35678a1b5316d212f4f366b2477232534a8aeca37f3c'
        def file = TestHelper.createInMemTempFile('foo', 'Hello world')
        expect:
        CacheHelper.hashFileSha256Impl0(file) == EXPECTED
        and:
        CacheHelper.hashFileSha256Impl0(file) == DigestUtils.sha256Hex(file.bytes)
    }

    def 'should hash content with sha256' () {
        given:
        def file = TestHelper.createInMemTempFile('foo.txt', 'Hello world')
        and:
        Global.session = Mock(Session) {
            getBaseDir() >> file.getParent()
            getCommitId() >> '123456'
        }

        expect:
        // this is not expecting to return the sha-256 of the file content BUT
        // the murmur hashing of the sha-256 string obtained by hashing the file
        CacheHelper.hasher(file, CacheHelper.HashMode.SHA256).hash().toString() == 'd29e7ba0fbcc617ab8e1e44e81381aed'
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
        def hash1 = CacheHelper.hashDirSha256(CacheHelper.defaultHasher().newHasher(), folder.resolve('dir1'), folder.resolve('dir1'))
        and:
        def hash2 = CacheHelper.hashDirSha256(CacheHelper.defaultHasher().newHasher(), folder.resolve('dir2'), folder.resolve('dir2'))

        then:
        hash1.hash() == hash2.hash()

    }
}
