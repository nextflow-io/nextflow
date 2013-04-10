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

import com.google.common.hash.Hashing
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CacheHelperTest extends Specification {

    def 'test hashCode' () {

        setup:
        def anInteger = Hashing.murmur3_32().newHasher().putInt(1).hash()
        def aLong = Hashing.murmur3_32().newHasher().putLong(100L).hash()
        def aString = Hashing.murmur3_32().newHasher().putString('Hola').hash()
        def aDouble = Hashing.murmur3_32().newHasher().putDouble(10.5).hash()
        def aBool = Hashing.murmur3_32().newHasher().putBoolean(true).hash()
        def aByteArray = Hashing.murmur3_32().newHasher().putBytes([0x1,0x2,0x3] as byte[]).hash()
        def anObjectArray = Hashing.murmur3_32().newHasher().putInt(1).putInt(2).putInt(3).hash()
        def aMap =  Hashing.murmur3_32().newHasher().putInt(1).putString('String1').putBoolean(true).hash()
        def aList = Hashing.murmur3_32().newHasher().putString('A').putString('B').putString('C').hash()

        def file = File.createTempFile('xxx','yyy')
        file.deleteOnExit()
        file.text = 'Hello'

        def aFile = Hashing.murmur3_32().newHasher().putString(file.absolutePath).putLong(file.length()).putLong(file.lastModified()).hash()


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
    }

    def 'test overflow' () {

        when:
        def hasher = Hashing.murmur3_32().newHasher().putInt(1)

        println hasher.hash()

        hasher = hasher.putInt(2)

        println hasher.hash()

        then:
        true

    }

}
