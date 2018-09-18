/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.util;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.Collection;
import java.util.Map;
import java.util.Set;
import java.util.UUID;

import com.google.common.hash.Funnels;
import com.google.common.hash.HashCode;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hasher;
import com.google.common.hash.Hashing;
import nextflow.extension.FilesEx;
import nextflow.file.FileHolder;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Provide helper method to handle caching
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class CacheHelper {

    public enum HashMode { STANDARD, DEEP }

    private static final Logger log = LoggerFactory.getLogger(CacheHelper.class);

    private static HashFunction DEFAULT_HASHING = Hashing.murmur3_128();

    private static int HASH_BITS = DEFAULT_HASHING.bits();

    private static int HASH_BYTES = HASH_BITS / 8;

    public static HashFunction defaultHasher() {
        return DEFAULT_HASHING;
    }

    public static Hasher hasher( Object value ) {
        return hasher(value, HashMode.STANDARD);
    }

    public static Hasher hasher( Object value, HashMode mode ) {
        return hasher( DEFAULT_HASHING, value, mode );
    }

    public static Hasher hasher( HashFunction function, Object value, HashMode mode ) {
        return hasher( function.newHasher(), value, mode );
    }

    public static Hasher hasher( Hasher hasher, Object value, HashMode mode ) {

        if( value == null )
            return hasher;

        if( value instanceof Boolean )
            return hasher.putBoolean((Boolean) value);

        if( value instanceof Short )
            return hasher.putShort((Short) value);

         if( value instanceof Integer)
            return hasher.putInt((Integer) value);

         if( value instanceof Long )
            return hasher.putLong((Long) value);

        if( value instanceof Float )
            return hasher.putFloat((Float) value);

        if( value instanceof Double )
            return hasher.putDouble( (Double)value );

        if( value instanceof Byte )
            return hasher.putByte( (Byte)value );

        if( value instanceof Number )
            // reduce all other number types (BigInteger, BigDecimal, AtomicXxx, etc) to string equivalent
            return hasher.putUnencodedChars(value.toString());

        if( value instanceof Character )
            return hasher.putChar( (Character)value );

        if( value instanceof CharSequence )
            return hasher.putUnencodedChars( (CharSequence)value );

        if( value instanceof byte[] )
            return hasher.putBytes( (byte[])value );

        if( value instanceof Object[]) {
            for( Object item: ((Object[])value) )
                hasher = CacheHelper.hasher( hasher, item, mode );
            return hasher;
        }

        if( value instanceof Map) {
            // note: should map be order invariant as Set ?
            for( Object item : ((Map)value).values() )
                hasher = CacheHelper.hasher( hasher, item, mode );
            return hasher;
        }

        if( value instanceof Bag || value instanceof Set )
            return hashUnorderedCollection(hasher, (Collection) value, mode);

        if( value instanceof Collection) {
            for( Object item: ((Collection)value) )
                hasher = CacheHelper.hasher( hasher, item, mode );
            return hasher;
        }

        if( value instanceof FileHolder )
            return CacheHelper.hasher(hasher, ((FileHolder) value).getSourceObj(), mode );

        if( value instanceof Path )
            return hashFile(hasher, (Path)value, mode);

        if( value instanceof java.io.File )
            return hashFile(hasher, (java.io.File)value, mode);

        if( value instanceof UUID ) {
            UUID uuid = (UUID)value;
            return hasher.putLong(uuid.getMostSignificantBits()).putLong(uuid.getLeastSignificantBits());
        }

        if( value instanceof VersionNumber ) {
            return hasher.putInt( value.hashCode() );
        }

        log.debug("[WARN] Unknown hashing type: {} -- {}", value.getClass(), value);
        return hasher.putInt( value.hashCode() );
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
    static private Hasher hashFile( Hasher hasher, java.io.File file, HashMode mode ) {
        return hashFile(hasher, file.toPath(), mode);
    }

    /**
     * Hashes the specified file
     *
     * @param hasher The current {@code Hasher} object
     * @param path The {@code Path} object to hash
     * @param mode When {@code mode} is equals to the string {@code deep} is used teh file content
     *   in order to create the hash key for this file, otherwise just the file metadata information
     *   (full name, size and last update timestamp)
     * @return The updated {@code Hasher} object
     */
    static private Hasher hashFile( Hasher hasher, Path path, HashMode mode ) {
        BasicFileAttributes attrs=null;
        try {
            attrs = Files.readAttributes(path, BasicFileAttributes.class);
        }
        catch(IOException e) {
            log.debug("Unable to get file attributes file: {} -- Cause: {}", path, e.toString());
        }

        if( mode==HashMode.DEEP && attrs!=null && attrs.isRegularFile() )
            return hashFileContent(hasher, path);
        else
            return hashFileMetadata(hasher, path, attrs);
    }

    /**
     * Hashes the file by using the metadata information: full path string, size and last update timestamp
     *
     * @param hasher The current {@code Hasher} object
     * @param file file The {@code Path} object to hash
     * @return The updated {@code Hasher} object
     */
    static private Hasher hashFileMetadata( Hasher hasher, Path file, BasicFileAttributes attrs ) {

        hasher = hasher.putUnencodedChars( file.toAbsolutePath().toString() );
        if( attrs != null ) {
            hasher = hasher.putLong(attrs.size());
            if( attrs.lastModifiedTime() != null) {
                hasher = hasher.putLong( attrs.lastModifiedTime().toMillis() );
            }
        }
        return hasher;
    }


    /**
     * Hashes the file by reading file content
     *
     * @param hasher The current {@code Hasher} object
     * @param path file The {@code Path} object to hash
     * @return The updated {@code Hasher} object
     */

    static private Hasher hashFileContent( Hasher hasher, Path path ) {

        OutputStream output = Funnels.asOutputStream(hasher);
        try {
            Files.copy(path, output);
        }
        catch( IOException e ) {
            throw new IllegalStateException("Unable to hash content: " + path, e);
        }
        finally {
            FilesEx.closeQuietly(output);
        }

        return hasher;
    }

    static HashCode hashContent( Path file ) {
        return hashContent(file, null);
    }

    static HashCode hashContent( Path file, HashFunction function ) {

        if( function == null )
            function = DEFAULT_HASHING;

        Hasher hasher = function.newHasher();
        return hashFileContent(hasher, file).hash();
    }

    static private Hasher hashUnorderedCollection(Hasher hasher, Collection collection, HashMode mode)  {

        byte[] resultBytes = new byte[HASH_BYTES];
        for (Object item : collection) {
            byte[] nextBytes = CacheHelper.hasher(item,mode).hash().asBytes();
            if( nextBytes.length != resultBytes.length )
                throw new IllegalStateException("All hash codes must have the same bit length");

            for (int i = 0; i < nextBytes.length; i++) {
                resultBytes[i] += nextBytes[i];
            }
        }

        return hasher.putBytes(resultBytes);

    }

}
