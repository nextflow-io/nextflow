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

import java.nio.file.Path

import com.google.common.hash.HashFunction
import com.google.common.hash.Hasher
import com.google.common.hash.Hashing
import groovy.util.logging.Slf4j
import nextflow.processor.FileHolder

/**
 * Provide helper method to handle caching
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class CacheHelper {

    static Hasher hasher( def value ) {
        hasher( Hashing.murmur3_32(), value )
    }

    static Hasher hasher( HashFunction function, def value ) {
        hasher( function.newHasher(), value )
    }

    static Hasher hasher( Hasher hasher, def value ) {
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
                    hasher = CacheHelper.hasher( hasher, item as Object )
                }
                break

            case Collection:
                value.each { hasher = CacheHelper.hasher( hasher, it )  }
                break

            case Map:
                value.each { name, item ->
                    hasher = CacheHelper.hasher( hasher, item )
                }
                break


            case FileHolder:
                hasher = CacheHelper.hasher(hasher, ((FileHolder) value).sourceObj )
                break

            case Path:
                hasher = hashPath(hasher, (Path)value)
                break;

            case File:
                hasher = hashFile(hasher, (File)value)
                break

            case UUID:
                def uuid = value as UUID
                hasher = hasher.putLong(uuid.mostSignificantBits).putLong(uuid.leastSignificantBits)
                break

            default:
                log.debug "Unknown hashing type: ${value.getClass()} -- ${value}"
                hasher = hasher.putInt( value.hashCode() )
        }


        return hasher
    }


    static private Hasher hashPath( Hasher hasher, Path file ) {

        hasher = hasher.putString( file.toAbsolutePath().normalize().toString() )

        if( file.exists() ) {
            return hasher.putLong( file.size() ) .putLong(file.lastModified())
        }
        else {
            return hasher
        }
    }

    static private Hasher hashFile( Hasher hasher, File file ) {

        hasher = hasher.putString(file.absolutePath)

        if( file.exists() ) {
            return hasher .putLong(file.length()).putLong( file.lastModified())
        }
        else {
            return hasher
        }


    }

}
