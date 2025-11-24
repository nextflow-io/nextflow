/*
 * Copyright 2013-2024, Seqera Labs
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

import static nextflow.Const.*;

import java.io.IOException;
import java.io.OutputStream;
import java.io.File;
import java.nio.file.FileVisitResult;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.ProviderMismatchException;
import java.nio.file.SimpleFileVisitor;
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
import nextflow.io.SerializableMarker;
import nextflow.script.types.Bag;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import static nextflow.util.CacheHelper.HashMode;


/**
 * Implements the hashing of objects
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
public class HashBuilder {

    private static final Logger log = LoggerFactory.getLogger(HashBuilder.class);

    private static final HashFunction DEFAULT_HASHING = Hashing.murmur3_128();

    private static final int HASH_BITS = DEFAULT_HASHING.bits();

    private static final int HASH_BYTES = HASH_BITS / 8;

    private static final Map<String,Object> FIRST_ONLY;

    static {
        FIRST_ONLY = new HashMap<>(1);
        FIRST_ONLY.put("firstOnly", Boolean.TRUE);
    }

    public static Hasher defaultHasher() {
        return DEFAULT_HASHING.newHasher();
    }

    private Hasher hasher = defaultHasher();

    private HashMode mode = HashMode.STANDARD;

    private Path basePath;

    public HashBuilder() {}

    public HashBuilder withHasher(Hasher hasher) {
        this.hasher = hasher;
        return this;
    }

    public HashBuilder withMode(HashMode mode) {
        this.mode = mode;
        return this;
    }

    public HashBuilder withBasePath(Path basePath) {
        this.basePath = basePath;
        return this;
    }

    public HashBuilder with(Object value) {

        if( value == null )
            return this;

        else if( value instanceof Boolean )
            hasher.putBoolean((Boolean) value);

        else if( value instanceof Short )
            hasher.putShort((Short) value);

        else if( value instanceof Integer)
            hasher.putInt((Integer) value);

        else if( value instanceof Long )
            hasher.putLong((Long) value);

        else if( value instanceof Float )
            hasher.putFloat((Float) value);

        else if( value instanceof Double )
            hasher.putDouble( (Double)value );

        else if( value instanceof Byte )
            hasher.putByte( (Byte)value );

        else if( value instanceof Number )
            // reduce all other number types (BigInteger, BigDecimal, AtomicXxx, etc) to string equivalent
            hasher.putUnencodedChars(value.toString());

        else if( value instanceof Character )
            hasher.putChar( (Character)value );

        else if( value instanceof CharSequence )
            hasher.putUnencodedChars( (CharSequence)value );

        else if( value instanceof byte[] )
            hasher.putBytes( (byte[])value );

        else if( value instanceof Object[])
            for( Object item : ((Object[])value) )
                with(item);

        // note: should map be order invariant as Set ?
        else if( value instanceof Map )
            for( Object item : ((Map)value).values() )
                with(item);

        else if( value instanceof Map.Entry ) {
            Map.Entry entry = (Map.Entry)value;
            with(entry.getKey());
            with(entry.getValue());
        }

        else if( value instanceof Bag || value instanceof Set )
            hashUnorderedCollection(hasher, (Collection) value, mode);

        else if( value instanceof Collection)
            for( Object item : ((Collection)value) )
                with(item);

        else if( value instanceof Path )
            hashFile(hasher, (Path)value, mode, basePath);

        else if( value instanceof java.io.File )
            hashFile(hasher, (java.io.File)value, mode, basePath);

        else if( value instanceof UUID ) {
            UUID uuid = (UUID)value;
            hasher.putLong(uuid.getMostSignificantBits()).putLong(uuid.getLeastSignificantBits());
        }

        else if( value instanceof VersionNumber )
            hasher.putInt( value.hashCode() );

        else if( value instanceof SerializableMarker)
            hasher.putInt( value.hashCode() );

        else if( value instanceof CacheFunnel )
            ((CacheFunnel)value).funnel(hasher, mode);

        else if( value instanceof Enum )
            hasher.putUnencodedChars( value.getClass().getName() + "." + value );

        else {
            Bolts.debug1(log, FIRST_ONLY, "[WARN] Unknown hashing type: " + value.getClass());
            hasher.putInt( value.hashCode() );
        }

        return this;
    }

    public Hasher getHasher() {
        return hasher;
    }

    public HashCode build() {
        return hasher.hash();
    }

    public static Hasher hasher( Hasher hasher, Object value, HashMode mode ) {

        return new HashBuilder()
                .withHasher(hasher)
                .withMode(mode)
                .with(value)
                .getHasher();
    }

    /**
     * Hash a file using only the relative file name instead of
     * the absolute file path.
     *
     * @param path
     * @param basePath
     * @param mode
     */
    public static HashCode hashPath(Path path, Path basePath, HashMode mode) {
        return new HashBuilder().withMode(mode).withBasePath(basePath).with(path).build();
    }

    /**
     * Hashes the specified file
     *
     * @param hasher The current {@code Hasher} object
     * @param file The {@code File} object to hash
     * @param mode When {@code mode} is equals to the string {@code deep} is used the file content
     *   in order to create the hash key for this file, otherwise just the file metadata information
     *   (full name, size and last update timestamp)
     * @return The updated {@code Hasher} object
     */
    static private Hasher hashFile( Hasher hasher, java.io.File file, HashMode mode, Path basePath ) {
        return hashFile(hasher, file.toPath(), mode, basePath);
    }

    /**
     * Hashes the specified file
     *
     * @param hasher The current {@code Hasher} object
     * @param path The {@code Path} object to hash
     * @param mode When {@code mode} is equals to the string {@code deep} is used the file content
     *   in order to create the hash key for this file, otherwise just the file metadata information
     *   (full name, size and last update timestamp)
     * @return The updated {@code Hasher} object
     */
    static private Hasher hashFile( Hasher hasher, Path path, HashMode mode, Path basePath ) {
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
        catch(Exception e) {
            log.warn("Unable to get file attributes file: {} -- Cause: {}", FilesEx.toUriString(path), e.toString());
        }

        if( (mode==HashMode.STANDARD || mode==HashMode.LENIENT) && isAssetFile(path, DEFAULT_ROOT) ) {
            if( attrs==null ) {
                // when file attributes are not avail, or it's a directory
                // hash the file using the file name path and the repository
                log.warn("Unable to fetch attribute for file: {} - Hash is inferred from Git repository commit Id", FilesEx.toUriString(path));
                return hashFileAsset(hasher, path);
            }
            final Path base = Global.getSession().getBaseDir();
            if( attrs.isDirectory() ) {
                // hash all the directory content
                return hashDirSha256(hasher, path, base);
            }
            else {
                // hash the content being an asset file
                // (i.e. included in the project repository) it's expected to small file
                // which makes the content hashing doable
                return hashFileSha256(hasher, path, base);
            }
        }

        if( mode==HashMode.DEEP && attrs!=null && attrs.isRegularFile() )
            return hashFileContent(hasher, path);
        if( mode==HashMode.SHA256 && attrs!=null && attrs.isRegularFile() )
            return hashFileSha256(hasher, path, null);
        // default
        return hashFileMetadata(hasher, path, attrs, mode, basePath);
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

    static protected Hasher hashFileSha256( Hasher hasher, Path path, Path base ) {
        try {
            log.trace("Hash sha-256 file content path={} - base={}", path, base);
            // the file relative base
            if( base!=null )
                hasher.putUnencodedChars(base.relativize(path).toString());
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

    /**
     * Compute an, order independent, hash of a directory path traversing recursively the directory content.
     *
     * @param hasher
     *      The {@link Hasher} object to which the resulting directory hash will be added.
     * @param dir
     *      The target directory path to be hashed.
     * @param base
     *      The "base" directory path against which resolve relative paths.
     * @return
     *      The resulting {@link Hasher} object updated with the directory path.
     */
    static protected Hasher hashDirSha256( Hasher hasher, Path dir, Path base ) {
        if( base==null )
            throw new IllegalArgumentException("Argument 'base' cannot be null");
        // the byte array used as "accumulator" for
        final byte[] resultBytes = new byte[HASH_BYTES];
        try {
            Files.walkFileTree(dir, new SimpleFileVisitor<Path>() {
                public FileVisitResult visitFile(Path path, BasicFileAttributes attrs) throws IOException {
                    log.trace("Hash sha-256 dir content [FILE] path={} - base={}", path, base);
                    try {
                        // the file relative base
                        final String relPath = base.relativize(path).toString();
                        // the file content sha-256 checksum
                        final String sha256 = sha256Cache.get(path);
                        // compute the file path hash and sum to the result hash
                        // since the sum is commutative, the traverse order does not matter
                        // compute a hash of the (file path, file hash) pair.
                        // since the sum is commutative, the resulting hash in `resultBytes` is invariant to the file traversal order.
                        // however, the file path and file hash do need to be processed together,
                        // otherwise this introduces an edge case with directories with similar contents with have the same sha (see nextflow-io/nextflow#6198)
                        sumBytes(resultBytes, hashBytes(Map.entry(relPath, sha256), HashMode.STANDARD));
                        return FileVisitResult.CONTINUE;
                    }
                    catch (ExecutionException t) {
                        throw new IOException(t);
                    }
                }

                public FileVisitResult preVisitDirectory(Path path, BasicFileAttributes attrs) {
                    log.trace("Hash sha-256 dir content [DIR] path={} - base={}", path, base);
                    // the file relative base
                    final String relPath = base.relativize(path).toString();
                    // compute the file path hash and sum to the result hash
                    // since the sum is commutative, the traverse order does not matter
                    sumBytes(resultBytes, hashBytes(relPath, HashMode.STANDARD));
                    return FileVisitResult.CONTINUE;
                }
            });
            // finally put the result bytes in the hashing
            hasher.putBytes(resultBytes);
        }
        catch (IOException t) {
            Throwable err = t.getCause()!=null ? t.getCause() : t;
            String msg = err.getMessage()!=null ? err.getMessage() : err.toString();
            log.warn("Unable to compute sha-256 hashing for directory: {} - Cause: {}", FilesEx.toUriString(dir), msg);
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
    static private Hasher hashFileMetadata( Hasher hasher, Path file, BasicFileAttributes attrs, HashMode mode, Path basePath ) {

        String filename = basePath != null && file.startsWith(basePath)
                ? basePath.relativize(file).toString()
                : file.toAbsolutePath().toString();

        hasher = hasher.putUnencodedChars( filename );
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
            // hash ghe collection item
            byte[] nextBytes = hashBytes(item, mode);
            // sum the hash bytes to the "resultBytes" accumulator
            // since the sum is a commutative operation the order does not matter
            sumBytes(resultBytes, nextBytes);
        }
        // add the result bytes and return the resulting object
        return hasher.putBytes(resultBytes);
    }

    static private byte[] hashBytes(Object item, HashMode mode) {
        return hasher(defaultHasher(), item, mode).hash().asBytes();
    }

    /**
     * Sum two arras of bytes having the same length, required to compute hash of unordered collections.
     *
     * - For each byte position, add the corresponding byte from nextBytes into resultBytes
     * - Order doesn't matter: addition is commutative (a + b = b + a), so the final result is
     *   the same no matter the order of items.
     * - This is what makes it suitable for unordered collections
     *
     * @param resultBytes
     *      The first argument to be summed. This array is used as the accumulator array (i.e. the result)
     * @param nextBytes
     *      The second argument to be summed.
     * @return
     *      The array resulting adding the bytes in the second array to the first one. Note,
     *      the result array instance is the same object passed as first argument.
     *
     */
    static private byte[] sumBytes(byte[] resultBytes, byte[] nextBytes) {
        if( nextBytes.length != resultBytes.length )
            throw new IllegalStateException("All hash codes must have the same bit length");
        for (int i = 0; i < nextBytes.length; i++) {
            resultBytes[i] += nextBytes[i];
        }
        return resultBytes;
    }

    /**
     * Check if the argument is an asset file i.e. a file that makes part of the
     * pipeline Git repository
     *
     * @param path
     *      The item to check.
     * @param assetRoot
     *      Location where assets are being stored.
     * @return
     *      Whether or not `path` is included in the pipeline Git repository.  
     */
    static protected boolean isAssetFile(Path path, File assetRoot) {
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
        return path.startsWith(session.getBaseDir()) || path.startsWith(assetRoot.toPath());
    }

}
