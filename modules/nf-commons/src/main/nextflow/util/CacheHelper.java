/*
 * Copyright 2020-2021, Seqera Labs
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

package nextflow.util;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.ProviderMismatchException;
import java.nio.file.attribute.BasicFileAttributes;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;
import java.util.UUID;
import java.util.concurrent.ExecutionException;

import com.google.common.cache.CacheBuilder;
import com.google.common.cache.CacheLoader;
import com.google.common.cache.LoadingCache;
import com.google.common.hash.Funnels;
import com.google.common.hash.HashCode;
import com.google.common.hash.HashFunction;
import com.google.common.hash.Hasher;
import com.google.common.hash.Hashing;
import com.google.common.io.ByteStreams;
import nextflow.Global;
import nextflow.ISession;
import nextflow.extension.Bolts;
import nextflow.extension.FilesEx;
import nextflow.file.FileHolder;
import nextflow.io.SerializableMarker;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * Provide helper method to handle caching
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class CacheHelper {

    public enum HashMode {

        STANDARD, DEEP, LENIENT, SHA256;

        public static HashMode of( Object obj ) {
            if( obj==null || obj instanceof Boolean )
                return null;
            if( obj instanceof CharSequence ) {
                if( "true".equals(obj) || "false".equals(obj) )
                    return null;
                if( "standard".equals(obj) )
                    return STANDARD;
                if( "lenient".equals(obj) )
                    return LENIENT;
                if( "deep".equals(obj) )
                    return DEEP;
                if( "sha256".equals(obj) )
                    return SHA256;
            }
            LoggerFactory.getLogger(HashMode.class).warn("Unknown cache mode: {}", obj.toString());
            return null;
        }
    }

    private static final Logger log = LoggerFactory.getLogger(CacheHelper.class);

    private static HashFunction DEFAULT_HASHING = Hashing.murmur3_128();

    private static int HASH_BITS = DEFAULT_HASHING.bits();

    private static int HASH_BYTES = HASH_BITS / 8;

    private static final Map<String,Object> FIRST_ONLY;

    static {
        FIRST_ONLY = new HashMap<>(1);
        FIRST_ONLY.put("firstOnly", Boolean.TRUE);
    }

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

        if( value instanceof Map ) {
            // note: should map be order invariant as Set ?
            for( Object item : ((Map)value).values() )
                hasher = CacheHelper.hasher( hasher, item, mode );
            return hasher;
        }

        if( value instanceof Map.Entry ) {
            Map.Entry entry = (Map.Entry)value;
            hasher = CacheHelper.hasher( hasher, entry.getKey(), mode );
            hasher = CacheHelper.hasher( hasher, entry.getValue(), mode );
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

        if( value instanceof SerializableMarker) {
            return hasher.putInt( value.hashCode() );
        }

        if( value instanceof CacheFunnel ) {
            return ((CacheFunnel) value).funnel(hasher,mode);
        }

        Bolts.debug1(log, FIRST_ONLY, "[WARN] Unknown hashing type: "+value.getClass());
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
            log.debug("Unable to get file attributes file: {} -- Cause: {}", FilesEx.toUriString(path), e.toString());
        }
        catch(ProviderMismatchException e) {
            // see https://github.com/nextflow-io/nextflow/pull/1382
            log.warn("File system is unable to get file attributes file: {} -- Cause: {}", FilesEx.toUriString(path), e.toString());
        }

        if( mode==HashMode.STANDARD && isAssetFile(path) ) {
            return attrs==null || attrs.isDirectory()
                    // when file attributes are not avail or it's a directory
                    // hash the file using the file name path and the repository
                    ? hashFileAsset(hasher, path)
                    // when it's not a directory, hash the content being an asset file
                    // (i.e. included in the project repository) it's expected to small file
                    // which makes the content hashing doable
                    : hashFileSha256(hasher, path);
        }

        if( mode==HashMode.DEEP && attrs!=null && attrs.isRegularFile() )
            return hashFileContent(hasher, path);
        if( mode==HashMode.SHA256 && attrs!=null && attrs.isRegularFile() )
            return hashFileSha256(hasher, path);
        // default
        return hashFileMetadata(hasher, path, attrs, mode);
    }


    static private LoadingCache<Path,String> sha256Cache = CacheBuilder
            .newBuilder()
            .maximumSize(10_000)
            .build(new CacheLoader<Path, String>() {
                @Override
                public String load(Path key) throws Exception {
                    return hashFileSha256Impl0(key);
                }
            });

    static private Hasher hashFileSha256( Hasher hasher, Path path ) {
        try {
            // file content hash
            String sha256 = sha256Cache.get(path);
            hasher.putUnencodedChars(sha256);
        }
        catch (ExecutionException t) {
            Throwable err = t.getCause()!=null ? t.getCause() : t;
            String msg = err.getMessage()!=null ? err.getMessage() : err.toString();
            log.warn("Unable to compute sha-256 hashing for file: {} - Cause: {}", FilesEx.toUriString(path), msg);
        }
        return hasher;
    }

    static protected String hashFileSha256Impl0(Path path) throws IOException {
        log.debug("Hash asset file sha-256: {}", path);
        Hasher hasher = Hashing.sha256().newHasher();
        ByteStreams.copy(Files.newInputStream(path), Funnels.asOutputStream(hasher));
        return hasher.hash().toString();
    }

    static private Hasher hashFileAsset( Hasher hasher, Path path ) {
        log.debug("Hash asset file: {}", path);
        hasher.putUnencodedChars( Global.getSession().getCommitId() );
        return hasher;
    }

    /**
     * Hashes the file by using the metadata information: full path string, size and last update timestamp
     *
     * @param hasher The current {@code Hasher} object
     * @param file file The {@code Path} object to hash
     * @return The updated {@code Hasher} object
     */
    static private Hasher hashFileMetadata( Hasher hasher, Path file, BasicFileAttributes attrs, HashMode mode ) {

        hasher = hasher.putUnencodedChars( file.toAbsolutePath().toString() );
        if( attrs != null ) {
            hasher = hasher.putLong(attrs.size());
            if( attrs.lastModifiedTime() != null && mode != HashMode.LENIENT ) {
                hasher = hasher.putLong( attrs.lastModifiedTime().toMillis() );
            }
        }

        if( log.isTraceEnabled() ) {
            log.trace("Hashing file meta: path={}; size={}, lastModified={}, mode={}",
                    file.toAbsolutePath().toString(),
                    attrs!=null ? attrs.size() : "--",
                    attrs!=null && attrs.lastModifiedTime() != null && mode != HashMode.LENIENT ? attrs.lastModifiedTime().toMillis() : "--",
                    mode
            );
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
            throw new IllegalStateException("Unable to hash content: " + FilesEx.toUriString(path), e);
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

    /**
     * Check if the argument is an asset file i.e. a file that makes part of the
     * pipeline Git repository
     *
     * @param path
     * @return
     */
    static protected boolean isAssetFile(Path path) {
        final ISession session = Global.getSession();
        if( session==null )
            return false;
        // if the commit ID is null the current run is not launched from a repo
        if( session.getCommitId()==null )
            return false;
        // if the file belong to different file system, cannot be a file belonging to the repo
        if( session.getBaseDir().getFileSystem()!=path.getFileSystem() )
            return false;
        // if the file is in the same directory as the base dir it's a asset by definition
        return path.startsWith(session.getBaseDir());
    }

}
