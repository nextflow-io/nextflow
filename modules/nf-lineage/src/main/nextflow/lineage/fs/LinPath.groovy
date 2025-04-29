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

    private String query

    private String fragment

    /*
     * Only needed to prevent serialization issues - see https://github.com/nextflow-io/nextflow/issues/5208
     */
    protected LinPath(){}

    LinPath(LinFileSystem fs, URI uri) {
        if( uri.scheme != SCHEME ) {
            throw new IllegalArgumentException("Invalid LID URI - scheme is different for $SCHEME")
        }
        this.fileSystem = fs
        setFieldsFormURI(uri)
        //Check if query and fragment are with filePath
        if (query == null && fragment == null){
            setFieldsFormURI(new URI(toUriString()))
        }
    }
    private void setFieldsFormURI(URI uri){
        this.query = uri.query
        this.fragment = uri.fragment
        this.filePath = resolve0(fileSystem, norm0("${uri.authority?:''}${uri.path}") )
    }

    protected LinPath(String query, String fragment, String filepath, LinFileSystem fs) {
        this.fileSystem = fs
        this.query = query
        this.fragment = fragment
        this.filePath = filepath
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

    protected static void validateDataOutput(FileOutput lidObject) {
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
            log.warn("Checksum of '$hashedPath' does not match with the one stored in the metadata")
    }

    protected static isAlgorithmSupported(String algorithm) {
        return algorithm && algorithm in SUPPORTED_CHECKSUM_ALGORITHMS
    }

    @TestOnly
    protected String getFilePath() { this.filePath }

    /**
     * Finds the target path of a LinPath.
     *
     * @param fs LinFileSystem associated to the LinPath to find
     * @param filePath Path associated to the LinPath to find
     * @param resultsAsPath True to return metadata descriptions as LinMetadataPath
     * @param children Sub-object/path inside the description
     * @return Path or LinMetadataPath associated to the LinPath
     * @throws Exception
     *      IllegalArgumentException if the filepath, filesystem or its LinStore are null.
     *      FileNotFoundException if the filePath or children are not found in the LinStore.
     */
    protected static Path findTarget(LinFileSystem fs, String filePath, boolean resultsAsPath, String[] children = []) throws Exception {
        if( !fs )
            throw new IllegalArgumentException("Cannot get target path for a relative lineage path")
        if( filePath.isEmpty() || filePath == SEPARATOR )
            throw new IllegalArgumentException("Cannot get target path for an empty lineage path")
        final store = fs.getStore()
        if( !store )
            throw new Exception("Lineage store not found - Check Nextflow configuration")
        final object = store.load(filePath)
        if( object ) {
            if( object instanceof FileOutput ) {
                return getTargetPathFromOutput(object, children)
            }
            if( resultsAsPath ) {
                return getMetadataAsTargetPath(object, fs, filePath, children)
            }
        } else {
            // If there isn't metadata check the parent to check if it is a subfolder of a task/workflow output
            final currentPath = Path.of(filePath)
            final parent = Path.of(filePath).getParent()
            if( parent ) {
                ArrayList<String> newChildren = new ArrayList<String>()
                newChildren.add(currentPath.getFileName().toString())
                newChildren.addAll(children)
                //resultsAsPath set to false because parent paths are only inspected for DataOutputs
                return findTarget(fs, parent.toString(), false, newChildren as String[])
            }
        }
        throw new FileNotFoundException("Target path '$filePath' does not exist")
    }

    protected static Path getMetadataAsTargetPath(LinSerializable results, LinFileSystem fs, String filePath, String[] children) {
        if( !results ) {
            throw new FileNotFoundException("Target path '$filePath' does not exist")
        }
        if( children && children.size() > 0 ) {
            return getSubObjectAsPath(fs, filePath, results, children)
        } else {
            return generateLinMetadataPath(fs, filePath, results, children)
        }
    }

    /**
     * Get a metadata sub-object as LinMetadataPath.
     * If the requested sub-object is the workflow or task outputs, retrieves the outputs from the outputs description.
     *
     * @param fs LinFilesystem for the te.
     * @param key Parent metadata key.
     * @param object Parent object.
     * @param children Array of string in indicating the properties to navigate to get the sub-object.
     * @return LinMetadataPath or null in it does not exist
     */
    static LinMetadataPath getSubObjectAsPath(LinFileSystem fs, String key, LinSerializable object, String[] children) {
        if( isSearchingOutputs(object, children) ) {
            // When asking for a Workflow or task output retrieve the outputs description
            final outputs = fs.store.load("${key}/output")
            if( !outputs ) {
                throw new FileNotFoundException("Target path '$key#output' does not exist")
            }
            return generateLinMetadataPath(fs, key, outputs, children)
        } else {
            return generateLinMetadataPath(fs, key, object, children)
        }
    }

    private static LinMetadataPath generateLinMetadataPath(LinFileSystem fs, String key, Object object, String[] children) {
        def creationTime = toFileTime(navigate(object, 'createdAt') as OffsetDateTime ?: OffsetDateTime.now())
        final output = children ? navigate(object, children.join('.')) : object
        if( !output ) {
            throw new FileNotFoundException("Target path '$key#${children.join('.')}' does not exist")
        }
        return new LinMetadataPath(encodeSearchOutputs(output, true), creationTime, fs, key, children)
    }

    private static Path getTargetPathFromOutput(FileOutput object, String[] children) {
        final lidObject = object as FileOutput
        // return the real path stored in the metadata
        validateDataOutput(lidObject)
        def realPath = FileHelper.toCanonicalPath(lidObject.path as String)
        if( children && children.size() > 0 )
            realPath = realPath.resolve(children.join(SEPARATOR))
        if( !realPath.exists() )
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
        return result ? new LinPath(fragment, query, result, null) : null
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
            return new LinPath(fragment, query, path.getName(index).toString(), null)
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
            return new LinPath(that.query, that.fragment, newPath.toString(), fileSystem)
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
        return new LinPath(lidOther.query, lidOther.fragment, path.getNameCount() > 0 ? path.toString() : SEPARATOR, null)
    }

    @Override
    URI toUri() {
        return asUri("${SCHEME}://${filePath}${query ? '?' + query : ''}${fragment ? '#' + fragment : ''}")
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
     * Get the path associated to a DataOutput metadata.
     *
     * @return Path associated to a DataOutput
     * @throws FileNotFoundException if the metadata associated to the LinPath does not exist or its type is not a DataOutput.
     */
    protected Path getTargetPath() {
        return findTarget(fileSystem, filePath, false, parseChildrenFromFragment(fragment))
    }

    /**
     * Get the path associated to any metadata object.
     *
     * @return Path associated to a DataOutput or LinMetadataFile with the metadata object for other types.
     * @throws FileNotFoundException if the metadata associated to the LinPath does not exist
     */
    protected Path getTargetOrMetadataPath() {
        return findTarget(fileSystem, filePath, true, parseChildrenFromFragment(fragment))
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
            throw new IllegalArgumentException("Invalid LID file system path URI - make sure the schema prefix does not container more than two slash characters - offending value: $path")
        if( path == LID_PROT ) //Empty path case
            return new URI("lid:///")
        return new URI(path)
    }

    @Override
    String toString() {
        return "$filePath${query ? '?' + query : ''}${fragment ? '#' + fragment : ''}".toString()
    }

}

