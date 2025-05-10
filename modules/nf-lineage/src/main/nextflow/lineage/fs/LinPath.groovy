/*
 * Copyright 2013-2025, Seqera Labs
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

package nextflow.lineage.fs

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.file.LogicalDataPath
import nextflow.lineage.model.Checksum
import nextflow.lineage.model.FileOutput
import nextflow.lineage.model.WorkflowRun
import nextflow.lineage.serde.LinSerializable
import nextflow.util.CacheHelper
import nextflow.util.TestOnly

import static LinFileSystemProvider.*
import static nextflow.lineage.LinUtils.*

import java.nio.file.FileSystem
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService
import java.time.OffsetDateTime

/**
 * LID file system path
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class LinPath implements Path, LogicalDataPath {

    static public final List<String> SUPPORTED_CHECKSUM_ALGORITHMS = ["nextflow"]
    static public final String SEPARATOR = '/'
    public static final String LID_PROT = "${SCHEME}://"

    static private final String[] EMPTY = new String[]{}

    private LinFileSystem fileSystem

    // String with the lineage file path
    private String filePath

    /*
     * Only needed to prevent serialization issues - see https://github.com/nextflow-io/nextflow/issues/5208
     */
    protected LinPath(){}

    LinPath(LinFileSystem fs, URI uri) {
        if( uri.scheme != SCHEME ) {
            throw new IllegalArgumentException("Invalid LID URI - scheme is different for $SCHEME")
        }
        this.fileSystem = fs
        this.filePath = resolve0(fileSystem, norm0("${uri.authority?:''}${uri.path}") )
        // Warn if query is specified
        if( uri.query )
            log.warn("Query string is not supported for Lineage URI: `$uri` -- it will be ignored")
        // Warn if fragment is specified
        if( uri.fragment )
            log.warn("Fragment is not supported for Lineage URI: `$uri` -- it will be ignored")
    }

    protected LinPath(String filePath, LinFileSystem fs) {
        this.fileSystem = fs
        this.filePath = filePath
    }

    LinPath(LinFileSystem fs, String path) {
        this( fs, asUri( LID_PROT + norm0(path)) )
    }

    LinPath(LinFileSystem fs, String first, String[] more) {
        this( fs, asUri( LID_PROT + buildPath(first, more) ) )
    }

    static String asUriString(String first, String... more) {
        return LID_PROT + buildPath(first, more)
    }

    static boolean isLidUri(String path) {
        return path && path.startsWith(LID_PROT)
    }

    private static String buildPath(String first, String[] more) {
        first = norm0(first)
        if( more ) {
            final morePath = norm0(more).join(SEPARATOR)
            return first.isEmpty() ? morePath : first + SEPARATOR + morePath
        }
        return first
    }

    protected static void validateFileOutput(FileOutput lidObject) {
        final hashedPath = FileHelper.toCanonicalPath(lidObject.path as String)
        if( !hashedPath.exists() )
            throw new FileNotFoundException("Target path $lidObject.path does not exist")
        validateChecksum(lidObject.checksum, hashedPath)
    }

    protected static void validateChecksum(Checksum checksum, Path hashedPath) {
        if( !checksum )
            return
        if( !isAlgorithmSupported(checksum.algorithm) ) {
            log.warn("Checksum of '$hashedPath' can't be validated. Algorithm '${checksum.algorithm}' is not supported")
            return
        }
        final hash = checksum.mode
            ? CacheHelper.hasher(hashedPath, CacheHelper.HashMode.of(checksum.mode.toString().toLowerCase())).hash().toString()
            : CacheHelper.hasher(hashedPath).hash().toString()
        if( hash != checksum.value )
            log.warn("Checksum of '$hashedPath' does not match with lineage metadata")
    }

    protected static isAlgorithmSupported(String algorithm) {
        return algorithm && algorithm in SUPPORTED_CHECKSUM_ALGORITHMS
    }

    @TestOnly
    protected String getFilePath() { this.filePath }

    protected List<Path> getSubPaths(){
        if( !fileSystem )
            throw new IllegalArgumentException("Cannot get sub-paths for a relative lineage path")
        if( filePath.isEmpty() || filePath == SEPARATOR )
            throw new IllegalArgumentException("Cannot get sub-paths for an empty lineage path (lid:///)")
        final store = fileSystem.getStore()
        if( !store )
            throw new Exception("Lineage store not found - Check Nextflow configuration")
        return store.getSubKeys(filePath).collect {new LinPath(fileSystem as LinFileSystem, it)} as List<Path>
    }

    /**
     * Finds the target path of a LinPath.
     * This method return the different types depending on the type of metadata pointing:
     * - When the LinPath point to FileOutput metadata or a subpath, it returns the real path.
     * - When it points other lineage metadata or a fragment of a lineage metadata and asMetadata is true, it returns a LinMetadataPath
     *   which contains in memory the lineage metadata description or the requested fragment of this description.
     * - When it points to a WorkflowRun or TaskRun metadata or subpath and asIntermediate is set to true. LinIntermediatePath, which is representing a directory
     * In other cases it will return a FileNotFoundException
     *
     * @param fs LinFileSystem associated to the LinPath to find
     * @param filePath Path to look for the target path
     * @param asMetadata Flag to indicate if other metadata descriptions must be returned as LinMetadataPath.
     * @param asIntermediate Flag to indicate if WorkflowRun and TaskRun subpaths must be returned as LinIntermediatePath.
     * @return Real Path, LinMetadataPath or LinIntermediatePath path associated to the LinPath
     * @throws Exception
     *      IllegalArgumentException if the filePath, filesystem or its LinStore are null.
     *      FileNotFoundException if the filePath is not found in the LinStore.
     */
    protected static Path findTarget(LinFileSystem fs, String filePath, boolean asMetadata, String subPath=null) throws Exception {
        if( !fs )
            throw new IllegalArgumentException("Cannot get target path for a relative lineage path")
        if( filePath.isEmpty() || filePath == SEPARATOR )
            throw new IllegalArgumentException("Cannot get target path for an empty lineage path (lid:///)")
        final store = fs.getStore()
        if( !store )
            throw new Exception("Lineage store not found - Check Nextflow configuration")
        final record = store.load(filePath)
        if( record instanceof FileOutput )
            return getFileOutputAsTargetPath(record, subPath)
        if( record && asMetadata )
            return getMetadataAsTargetPath(record, fs, filePath)
        // recursively check parent paths for metadata descriptions
        final currentPath = Path.of(filePath)
        final parent = currentPath.getParent()
        if( parent ) {
            final filename = currentPath.getFileName().toString()
            subPath = subPath
                ? "${filename}${SEPARATOR}${subPath}".toString()
                : filename
            return findTarget(fs, parent.toString(), false, subPath)
        }
        throw new FileNotFoundException("Target path '${filePath}' does not exist")
    }

    protected static Path getMetadataAsTargetPath(record, LinFileSystem fs, String filePath) {
        if( !record )
            throw new FileNotFoundException("Target path '$filePath' does not exist")
        final creationTime = toFileTime(navigate(record, 'createdAt') as OffsetDateTime ?: OffsetDateTime.now())
        return new LinMetadataPath(encodeSearchOutputs(record, true), creationTime, fs, filePath)
    }

    private static Path getFileOutputAsTargetPath(FileOutput record, String subPath) {
        // return the real path stored in the metadata
        validateFileOutput(record)
        def realPath = FileHelper.toCanonicalPath(record.path as String)
        if( subPath )
            realPath = realPath.resolve(subPath)
        if (!realPath.exists())
            throw new FileNotFoundException("Target path '$realPath' does not exist")
        return realPath
    }

    private static boolean isEmptyBase(LinFileSystem fs, String base) {
        return !base || base == SEPARATOR || (fs && base == "..")
    }

    private static String resolve0(LinFileSystem fs, String base, String[] more) {
        if( isEmptyBase(fs, base) ) {
            return resolveEmptyPathCase(fs, more as List)
        }
        if( base.contains(SEPARATOR) ) {
            final parts = base.tokenize(SEPARATOR)
            final remain = parts[1..-1] + more.toList()
            return resolve0(fs, parts[0], remain as String[])
        }
        final result = Path.of(base)
        return more ? result.resolve(more.join(SEPARATOR)).toString() : result.toString()
    }

    private static String resolveEmptyPathCase(LinFileSystem fs, List<String> more) {
        switch( more.size() ) {
            case 0:
                return "/"
            case 1:
                return resolve0(fs, more[0], EMPTY)
            default:
                return resolve0(fs, more[0], more[1..-1] as String[])
        }
    }

    static private String norm0(String path) {
        if( !path || path == SEPARATOR )
            return ""
        //Remove repeated elements
        path = Path.of(path.trim()).normalize().toString()
        //Remove initial and final separators
        if( path.startsWith(SEPARATOR) )
            path = path.substring(1)
        if( path.endsWith(SEPARATOR) )
            path = path.substring(0, path.size() - 1)
        return path
    }

    static private String[] norm0(String... path) {
        for( int i = 0; i < path.length; i++ ) {
            path[i] = norm0(path[i])
        }
        return path
    }

    @Override
    FileSystem getFileSystem() {
        return fileSystem
    }

    @Override
    boolean isAbsolute() {
        return fileSystem != null
    }

    @Override
    Path getRoot() {
        return new LinPath(fileSystem, SEPARATOR)
    }

    @Override
    Path getFileName() {
        final result = Path.of(filePath).getFileName()?.toString()
        return result ? new LinPath(result, null) : null
    }

    @Override
    Path getParent() {
        final c = getNameCount()
        if( c > 1 )
            return subpath(0, c - 1)
        if( c == 1 )
            return new LinPath(fileSystem, SEPARATOR)
        return null
    }

    @Override
    int getNameCount() {
        return Path.of(filePath).nameCount
    }

    @Override
    Path getName(int index) {
        if( index < 0 )
            throw new IllegalArgumentException("Path name index cannot be less than zero - offending value: $index")
        final path = Path.of(filePath)
        if( index == path.nameCount - 1 ) {
            return new LinPath(path.getName(index).toString(), null)
        }
        return new LinPath(index == 0 ? fileSystem : null, path.getName(index).toString())
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        if( beginIndex < 0 )
            throw new IllegalArgumentException("subpath begin index cannot be less than zero - offending value: $beginIndex")
        final path = Path.of(filePath)
        return new LinPath(beginIndex == 0 ? fileSystem : null, path.subpath(beginIndex, endIndex).toString())
    }

    @Override
    Path normalize() {
        return new LinPath(fileSystem, Path.of(filePath).normalize().toString())
    }

    @Override
    boolean startsWith(Path other) {
        return startsWith(other.toString())
    }

    @Override
    boolean startsWith(String other) {
        return filePath.startsWith(other)
    }

    @Override
    boolean endsWith(Path other) {
        return endsWith(other.toString())
    }

    @Override
    boolean endsWith(String other) {
        return filePath.endsWith(other)
    }

    @Override
    Path resolve(Path other) {
        if( LinPath.class != other.class )
            throw new ProviderMismatchException()

        final that = (LinPath) other

        if( that.fileSystem && this.fileSystem != that.fileSystem )
            return other
        if( that.isAbsolute() ) {
            return that
        } else {
            final newPath = Path.of(filePath).resolve(that.toString())
            return new LinPath(newPath.toString(), fileSystem)
        }
    }

    @Override
    Path resolve(String path) {
        if( !path )
            return this
        final scheme = FileHelper.getUrlProtocol(path)
        if( !scheme ) {
            // consider the path as a lid relative path
            return resolve(new LinPath(null, path))
        }
        if( scheme != SCHEME ) {
            throw new ProviderMismatchException()
        }
        final that = fileSystem.provider().getPath(asUri(path))
        return resolve(that)
    }

    @Override
    Path relativize(Path other) {
        if( LinPath.class != other.class ) {
            throw new ProviderMismatchException()
        }
        LinPath lidOther = other as LinPath
        if( this.isAbsolute() != lidOther.isAbsolute() )
            throw new IllegalArgumentException("Cannot compare absolute with relative paths");
        def path
        if( this.isAbsolute() ) {
            // Compare 'filePath' as absolute paths adding the root separator
            path = Path.of(SEPARATOR + filePath).relativize(Path.of(SEPARATOR + lidOther.filePath))
        } else {
            // Compare 'filePath' as relative paths
            path = Path.of(filePath).relativize(Path.of(lidOther.filePath))
        }
        return new LinPath(path.getNameCount() > 0 ? path.toString() : SEPARATOR, null)
    }

    @Override
    URI toUri() {
        return asUri("${SCHEME}://${filePath}")
    }

    String toUriString() {
        return toUri().toString()
    }

    @Override
    Path toAbsolutePath() {
        return this
    }

    @Override
    Path toRealPath(LinkOption... options) throws IOException {
        return this.getTargetOrMetadataPath()
    }

    Path toTargetPath() {
        return getTargetOrMetadataPath()
    }

    /**
     * Get the path associated with a FileOutput record.
     *
     * @return Path associated with a FileOutput record
     * @throws FileNotFoundException if the record does not exist or its type is not a FileOutput.
     */
    protected Path getTargetPath() {
        return findTarget(fileSystem, filePath, false, false)
    }

    /**
     * Get the path associated with a FileOutput record or an intermediate subpath.
     *
     * @return Path associated with a FileOutput record or a LinIntermediatePath if LinPath points to a workflow and task run subpath.
     * @throws FileNotFoundException if the record does not exist or its type is not a FileOutput or a intermediate directory
     */
    protected Path getTargetOrIntermediatePath() {
        return findTarget(fileSystem, filePath, false, true)
    }

    /**
     * Get the path associated with a lineage record.
     *
     * @return Path associated with a FileOutput record or a LinMetadataFile with the lineage record for other types, or a intermediate directory
     * @throws FileNotFoundException if the record does not exist
     */
    protected Path getTargetOrMetadataPath() {
        return findTarget(fileSystem, filePath, true, false)
    }

    @Override
    File toFile() throws IOException {
        throw new UnsupportedOperationException("toFile not supported by LinPath")
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) throws IOException {
        throw new UnsupportedOperationException("Register not supported by LinPath")
    }

    @Override
    int compareTo(Path other) {
        return toString().compareTo(other.toString());
    }

    @Override
    boolean equals(Object other) {
        if( LinPath.class != other.class ) {
            return false
        }
        final that = (LinPath) other
        return this.fileSystem == that.fileSystem && this.filePath.equals(that.filePath)
    }

    /**
     * @return The unique hash code for this path
     */
    @Override
    int hashCode() {
        return Objects.hash(fileSystem, filePath)
    }

    static URI asUri(String path) {
        if( !path )
            throw new IllegalArgumentException("Missing 'path' argument")
        if( !path.startsWith(LID_PROT) )
            throw new IllegalArgumentException("Invalid LID file system path URI - it must start with '${LID_PROT}' prefix - offendinf value: $path")
        if( path.startsWith(LID_PROT + SEPARATOR) && path.length() > 7 )
            throw new IllegalArgumentException("Invalid LID file system path URI - make sure the schema prefix does not container more than two slash characters or a query in the root '/' - offending value: $path")
        if( path == LID_PROT ) //Empty path case
            return new URI("lid:///")
        return new URI(path)
    }

    @Override
    String toString() {
        return filePath
    }

}

