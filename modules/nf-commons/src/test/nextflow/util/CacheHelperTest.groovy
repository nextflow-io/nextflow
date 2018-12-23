/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

import spock.lang.Specification

import java.nio.file.Files
import java.nio.file.attribute.FileTime

import com.google.common.hash.Hashing
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
        VALUE       | EXPECTED
        null        | null
        true        | null
        'true'      | null
        'false'     | null
        'standard'  | CacheHelper.HashMode.STANDARD
        'deep'      | CacheHelper.HashMode.DEEP
        'lenient'   | CacheHelper.HashMode.LENIENT

    }
}
