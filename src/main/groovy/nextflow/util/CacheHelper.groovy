/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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
import java.nio.file.Path

import embed.com.google.common.hash.Funnels
import embed.com.google.common.hash.HashCode
import embed.com.google.common.hash.HashFunction
import embed.com.google.common.hash.Hasher
import embed.com.google.common.hash.Hashing
import groovy.util.logging.Slf4j
import nextflow.processor.FileHolder

enum HashMode { STANDARD, DEEP }

/**
 * Provide helper method to handle caching
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class CacheHelper {

    static HashFunction defaultHasher() {
        Hashing.murmur3_128()
    }

    static Hasher hasher( value, HashMode mode = null) {
        hasher( defaultHasher(), value, mode )
    }

    static Hasher hasher( HashFunction function, value, HashMode mode = null ) {
        hasher( function.newHasher(), value, mode )
    }

    static Hasher hasher( Hasher hasher, value, HashMode mode = null ) {
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
                hasher = hasher.putUnencodedChars( value as CharSequence );
                break

            case Byte:
                hasher = hasher.putByte( value as Byte );
                break

            case (byte[]):
                hasher = hasher.putBytes( value as byte[] );
                break

            case (Object[]):
                for( def item: (value as Object[]) ) {
                    hasher = CacheHelper.hasher( hasher, item as Object, mode )
                }
                break

            case Collection:
                value.each { hasher = CacheHelper.hasher( hasher, it, mode )  }
                break

            case Map:
                value.each { name, item ->
                    hasher = CacheHelper.hasher( hasher, item, mode )
                }
                break


            case FileHolder:
                hasher = CacheHelper.hasher(hasher, ((FileHolder) value).sourceObj, mode )
                break

            case Path:
                hasher = hashFile(hasher, (Path)value, mode)
                break;

            case File:
                hasher = hashFile(hasher, (File)value, mode)
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

    /**
     * Hashes the specified file
     *
     * @param hasher The current {@code Hasher} object
     * @param file The {@code File} object to hash
     * @param mode When {@code mode} is equals to the string {@code deep} is used teh file content
     *   in order to create the hash key for this file, otherwise just the file metadata information
     *   (full name, size and last update timestamp)
     * @return The updated {@code Hasher} object
     */
    static private Hasher hashFile( Hasher hasher, File file, HashMode mode = null ) {
        hashFile(hasher, file.toPath(), mode)
    }

    /**
     * Hashes the specified file
     *
     * @param hasher The current {@code Hasher} object
     * @param file The {@code Path} object to hash
     * @param mode When {@code mode} is equals to the string {@code deep} is used teh file content
     *   in order to create the hash key for this file, otherwise just the file metadata information
     *   (full name, size and last update timestamp)
     * @return The updated {@code Hasher} object
     */
    static private Hasher hashFile( Hasher hasher, Path path, HashMode mode = null ) {
        if( mode == HashMode.DEEP && Files.isRegularFile(path))
            hashFileContent(hasher, path)

        else
            hashFileMetadata(hasher, path)
    }

    /**
     * Hashes the file by using the metadata information: full path string, size and last update timestamp
     *
     * @param hasher The current {@code Hasher} object
     * @param file file The {@code Path} object to hash
     * @return The updated {@code Hasher} object
     */
    static private Hasher hashFileMetadata( Hasher hasher, Path file ) {

        hasher = hasher.putUnencodedChars( file.toAbsolutePath().normalize().toString() )

        if( file.exists() ) {
            return hasher.putLong( file.size() ) .putLong(file.lastModified())
        }
        else {
            return hasher
        }
    }


    /**
     * Hashes the file by reading file content
     *
     * @param hasher The current {@code Hasher} object
     * @param file file The {@code Path} object to hash
     * @return The updated {@code Hasher} object
     */

    static private Hasher hashFileContent( Hasher hasher, Path path ) {

        def output = Funnels.asOutputStream(hasher)
        path.withInputStream { input ->
            output << input
        }
        output.closeQuietly()
        hasher
    }


    static HashCode hashContent( Path file, HashFunction function = null ) {

        if( !function )
            function = defaultHasher()

        Hasher hasher = function.newHasher();
        hashFileContent(hasher, file).hash()
    }

}
