/*
 * Copyright (c) 2012, the authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.util

import java.nio.file.Files

import embed.com.google.common.hash.Hashing
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CacheHelperTest extends Specification {

    def 'test hashCode' () {

        setup:
        def anInteger = Hashing.murmur3_128().newHasher().putInt(1).hash()
        def aLong = Hashing.murmur3_128().newHasher().putLong(100L).hash()
        def aString = Hashing.murmur3_128().newHasher().putString('Hola').hash()
        def aDouble = Hashing.murmur3_128().newHasher().putDouble(10.5).hash()
        def aBool = Hashing.murmur3_128().newHasher().putBoolean(true).hash()
        def aByteArray = Hashing.murmur3_128().newHasher().putBytes([0x1,0x2,0x3] as byte[]).hash()
        def anObjectArray = Hashing.murmur3_128().newHasher().putInt(1).putInt(2).putInt(3).hash()
        def aMap =  Hashing.murmur3_128().newHasher().putInt(1).putString('String1').putBoolean(true).hash()
        def aList = Hashing.murmur3_128().newHasher().putString('A').putString('B').putString('C').hash()

        def file = File.createTempFile('xxx','yyy')
        file.deleteOnExit()
        file.text = 'Hello'

        def aFile = Hashing.murmur3_128().newHasher().putString(file.absolutePath).putLong(file.length()).putLong(file.lastModified()).hash()


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



}
