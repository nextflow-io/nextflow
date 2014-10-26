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

import com.google.common.hash.Funnels
import com.google.common.hash.HashCode
import com.google.common.hash.HashFunction
import com.google.common.hash.Hasher
import com.google.common.hash.Hashing
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.extension.FilesEx
import nextflow.file.FileHolder

enum HashMode { STANDARD, DEEP }

/**
 * Provide helper method to handle caching
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
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

        /*
         * Used huge IF instead of Switch because the groovy switch implementation is very slow
         */

        if( value == null )
            return hasher

        if( value instanceof Boolean )
            return hasher.putBoolean( (Boolean)value )

        if( value instanceof Short )
            return hasher.putShort( (Short)value )

         if( value instanceof Integer)
            return hasher.putInt( (Integer)value );

         if( value instanceof Long )
            return hasher.putLong( (Long)value )

        if( value instanceof Float )
            return hasher.putFloat( (Float)value );

        if( value instanceof Double )
            return hasher.putDouble( (Double)value );

        if( value instanceof Character )
            return hasher.putChar( (Character)value );

        if( value instanceof CharSequence )
            return hasher.putUnencodedChars( (CharSequence)value );

        if( value instanceof Byte )
            return hasher.putByte( (Byte)value );

        if( value instanceof byte[] )
            return hasher.putBytes( (byte[])value );

        if( value instanceof Object[]) {
            for( Object item: ((Object[])value) )
                hasher = CacheHelper.hasher( hasher, item as Object, mode )
            return hasher
        }

        if( value instanceof Collection ) {
            for( Object item: ((Collection)value) )
                hasher = CacheHelper.hasher( hasher, item as Object, mode )
            return hasher
        }

        if( value instanceof Map ) {
            for( Object item : ((Map)value).values() )
                hasher = CacheHelper.hasher( hasher, item as Object, mode )
            return hasher
        }

        if( value instanceof FileHolder )
            return  CacheHelper.hasher(hasher, ((FileHolder) value).sourceObj, mode )

        if( value instanceof Path )
            return hashFile(hasher, (Path)value, mode)

        if( value instanceof  File )
            return hashFile(hasher, (File)value, mode)

        if( value instanceof UUID ) {
            def uuid = value as UUID
            return hasher.putLong(uuid.mostSignificantBits).putLong(uuid.leastSignificantBits)
        }

        log.debug "Unknown hashing type: ${value.getClass()} -- ${value}"
        return hasher.putInt( value.hashCode() )
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

        if( FilesEx.exists(file) ) {
            return hasher.putLong(file.size()) .putLong(FilesEx.lastModified(file))
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
        FilesEx.closeQuietly(output)
        hasher
    }


    static HashCode hashContent( Path file, HashFunction function = null ) {

        if( !function )
            function = defaultHasher()

        Hasher hasher = function.newHasher();
        hashFileContent(hasher, file).hash()
    }

}
