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

import com.google.common.hash.HashFunction
import com.google.common.hash.Hasher
import com.google.common.hash.Hashing
import groovy.util.logging.Slf4j
/**
 * Provide helper method to handle caching
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class CacheHelper {

    static Hasher hasher( def value ) {
        hasher( value, Hashing.murmur3_32() )
    }

    static Hasher hasher( def value, HashFunction function ) {
        hasher( value, function.newHasher() )
    }

    static Hasher hasher( def value, Hasher hasher ) {
        assert hasher

        if( value == null ) return hasher

        switch (value.getClass()) {
            case Boolean:
                hasher = hasher.putBoolean( value as Boolean );
                break

            case Short:
                hasher = hasher.putShort( value as Short );
                break

            case Integer:
                hasher = hasher.putInt( value as Integer );
                break

            case Long:
                hasher = hasher.putLong( value as Long );
                break

            case Float:
                hasher = hasher.putFloat( value as Float );
                break

            case Double:
                hasher = hasher.putDouble( value as Double );
                break

            case Character:
                hasher = hasher.putChar( value as Character );
                break

            case CharSequence:
                hasher = hasher.putString( value as CharSequence );
                break

            case Byte:
                hasher = hasher.putByte( value as Byte );
                break

            case (byte[]):
                hasher = hasher.putBytes( value as byte[] );
                break

            case (Object[]):
                for( def item: (value as Object[]) ) {
                    hasher = CacheHelper.hasher( item as Object, hasher )
                }
                break

            case Collection:
                value.each { hasher = CacheHelper.hasher( it,hasher )  }
                break

            case Map:
                value.each { name, item ->
                    hasher = CacheHelper.hasher( item, hasher )
                }
                break

            case File:
                def file = value as File
                hasher
                    .putString( file.absolutePath )
                    .putLong( file.length() )
                    .putLong( file.lastModified() )
                break

            default:
                log.debug "Unknown hashing type: ${value.getClass()} -- ${value}"
                hasher = hasher.putInt( value.hashCode() )
        }


        return hasher
    }



}
