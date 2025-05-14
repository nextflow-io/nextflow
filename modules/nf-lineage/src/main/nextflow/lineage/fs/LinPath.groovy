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

import static nextflow.lineage.LinUtils.*
import static nextflow.lineage.fs.LinFileSystemProvider.*

import java.nio.file.FileSystem
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService
import java.time.OffsetDateTime
import java.util.stream.Stream

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.file.FileHelper
import nextflow.file.LogicalDataPath
import nextflow.lineage.LinPropertyValidator
import nextflow.lineage.LinStore
import nextflow.lineage.model.v1beta1.Checksum
import nextflow.lineage.model.v1beta1.FileOutput
import nextflow.lineage.model.v1beta1.TaskRun
import nextflow.lineage.model.v1beta1.WorkflowRun
import nextflow.lineage.serde.LinSerializable
import nextflow.util.CacheHelper
import nextflow.util.TestOnly
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
        // Check if query and fragment are with filePath
        if( query == null && fragment == null )
            setFieldsFormURI(new URI(toUriString()))
        // Warn if query is specified
        if( query )
            log.warn("Query string is not supported for Lineage URI: `$uri` -- it will be ignored")
        // Validate fragment
        if( fragment )
            new LinPropertyValidator().validate(fragment.tokenize('.'))
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

    @Memoized
    protected static String validateDataOutput(FileOutput lidObject) {
        final hashedPath = FileHelper.toCanonicalPath(lidObject.path as String)
        if( !hashedPath.exists() )
            throw new FileNotFoundException("Target path $lidObject.path does not exist")
        return validateChecksum(lidObject.checksum, hashedPath)
    }

    protected static String validateChecksum(Checksum checksum, Path hashedPath) {
        if( !checksum )
            return null
        if( !isAlgorithmSupported(checksum.algorithm) ) {
            return "Checksum of '$hashedPath' can't be validated - algorithm '${checksum.algorithm}' is not supported"
        }
        final hash = checksum.mode
            ? CacheHelper.hasher(hashedPath, CacheHelper.HashMode.of(checksum.mode.toString().toLowerCase())).hash().toString()
            : CacheHelper.hasher(hashedPath).hash().toString()
        return hash != checksum.value
            ? "Checksum of '$hashedPath' does not match with lineage metadata"
            : null
    }

    protected static isAlgorithmSupported(String algorithm) {
        return algorithm && algorithm in SUPPORTED_CHECKSUM_ALGORITHMS
    }

    @TestOnly
    protected String getFilePath() { this.filePath }

    protected Stream<Path> getSubPaths(){
        if( !fileSystem )
            throw new IllegalArgumentException("Cannot get sub-paths for a relative lineage path")
        if( filePath.isEmpty() || filePath == SEPARATOR )
            throw new IllegalArgumentException("Cannot get sub-paths for an empty lineage path (lid:///)")
        final store = fileSystem.getStore()
        if( !store )
            throw new Exception("Lineage store not found - Check Nextflow configuration")
        return store.getSubKeys(filePath).map {new LinPath(fileSystem as LinFileSystem, it) as Path }
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
     * @param fragment String with path to sub-object inside the description
     * @param asMetadata Flag to indicate if other metadata descriptions must be returned as LinMetadataPath.
     * @param asIntermediate Flag to indicate if WorkflowRun and TaskRun subpaths must be returned as LinIntermediatePath.
     * @param subpath subpath associated to the target path to find. Used when looking for a parent path
     * @return Real Path, LinMetadataPath or LinIntermediatePath path associated to the LinPath
     * @throws Exception
     *      IllegalArgumentException if the filepath, filesystem or its LinStore are null.
     *      FileNotFoundException if the filePath, subpath and fragment is not found.
     */
    protected static Path findTarget(LinFileSystem fs, String filePath, String fragment, boolean asMetadata, boolean asIntermediate) throws Exception {
        if( !fs )
            throw new IllegalArgumentException("Cannot get target path for a relative lineage path")
        if( filePath.isEmpty() || filePath == SEPARATOR )
            throw new IllegalArgumentException("Cannot get target path for an empty lineage path (lid:///)")
        final store = fs.getStore()
        if( !store )
            throw new Exception("Lineage store not found - Check Nextflow configuration")
        findTarget0(fs, store, filePath, fragment, asMetadata, asIntermediate, [])
    }
    
    private static Path findTarget0(LinFileSystem fs, LinStore store, String filePath, String fragment, boolean asMetadata, boolean asIntermediate, List<String> subpath) {
        final object = store.load(filePath)
        if( object ) {
            return getTargetPathFromObject(object, fs, filePath, fragment, asMetadata, asIntermediate, subpath)
        } else {
            if( fragment ) {
                // If object doesn't exit, it's not possible to get fragment.
                throw new FileNotFoundException("Target path '$filePath#$fragment' does not exist")
            }
            return findTargetFromParent(fs, store, filePath, asIntermediate, subpath)
        }
    }

    private static Path findTargetFromParent(LinFileSystem fs, LinStore store, String filePath, boolean asIntermediate, List<String> subpath) {
        final currentPath = Path.of(filePath)
        final parent = Path.of(filePath).getParent()
        if( !parent ) {
            throw new FileNotFoundException("Target path '$filePath/${subpath.join('/')} does not exist")
        }
        ArrayList<String> newChildren = new ArrayList<String>()
        newChildren.add(currentPath.getFileName().toString())
        newChildren.addAll(subpath)
        //As Metadata set as false because parent path only inspected for FileOutput or intermediate.
        return findTarget0(fs, store, parent.toString(), null, false, asIntermediate, newChildren)
    }

    private static Path getTargetPathFromObject(LinSerializable object, LinFileSystem fs, String filePath, String fragment, boolean asMetadataPath, boolean asIntermediatePath,List<String> subpath) {
        // It's not possible to get a target path with both fragment and subpath
        if( fragment && subpath ) {
            throw new FileNotFoundException("Unable to get a target path for '$filePath' with fragments and subpath")
        }
        // If metadata flag is active and looks for a fragment returns the metadata despite the type of object
        if( asMetadataPath && fragment ){
           return getMetadataAsTargetPath(object, fs, filePath, fragment)
        }
        // Return real files when FileOutput sub-path
        if( object instanceof FileOutput ) {
            return getTargetPathFromOutput(object, subpath)
        }
        // Intermediate run case
        if( asIntermediatePath && (object instanceof WorkflowRun || object instanceof TaskRun) ) {
            return new LinIntermediatePath(fs, "$filePath/${subpath.join('/')}")
        }

        // It is not possible to get a metadata path with subpath. For other cases return metadata path if activated or throw exception
        if( asMetadataPath && !subpath)
            return getMetadataAsTargetPath(object, fs, filePath, fragment)
        else
            throw new FileNotFoundException("Target path '${filePath}/${subpath ? '/' + subpath.join('/') : ''}${fragment ? '#' + fragment : ''}' does not exist")
    }

    protected static Path getMetadataAsTargetPath(LinSerializable results, LinFileSystem fs, String filePath, String fragment) {
        if( !results ) {
            throw new FileNotFoundException("Target path '$filePath' does not exist")
        }
        if( fragment ) {
            return getSubObjectAsPath(fs, filePath, results, fragment)
        } else {
            return generateLinMetadataPath(fs, filePath, results, fragment)
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
    static LinMetadataPath getSubObjectAsPath(LinFileSystem fs, String key, LinSerializable object, String fragment) {
        if( isSearchingOutputs(object, fragment) ) {
            // When asking for a Workflow or task output retrieve the outputs description
            final outputs = fs.store.load("${key}#output")
            if( !outputs ) {
                throw new FileNotFoundException("Target path '$key#output' does not exist")
            }
            return generateLinMetadataPath(fs, key, outputs, fragment)
        } else {
            return generateLinMetadataPath(fs, key, object, fragment)
        }
    }

    private static LinMetadataPath generateLinMetadataPath(LinFileSystem fs, String key, Object object, String fragment) {
        def creationTime = toFileTime(navigate(object, 'createdAt') as OffsetDateTime ?: OffsetDateTime.now())
        final output = fragment ? navigate(object, fragment) : object
        if( !output ) {
            throw new FileNotFoundException("Target path '$key#${fragment}' does not exist")
        }
        return new LinMetadataPath(encodeSearchOutputs(output, true), creationTime, fs, key, fragment)
    }

    private static Path getTargetPathFromOutput(FileOutput object, List<String> children) {
        final lidObject = object as FileOutput
        // verify checksum validation
        final violation = validateDataOutput(lidObject)
        if( violation )
            log.warn1(violation)
        // return the real path stored in the metadata
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
        return result ? new LinPath(query, fragment, result, null) : null
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
            return new LinPath( query, fragment, path.getName(index).toString(), null)
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
     * Get the path associated with a FileOutput record.
     *
     * @return Path associated with a FileOutput record
     * @throws FileNotFoundException if the record does not exist or its type is not a FileOutput.
     */
    protected Path getTargetPath() {
        return findTarget(fileSystem, filePath, fragment, false, false)
    }

    /**
     * Get the path associated with a FileOutput record or an intermediate subpath.
     *
     * @return Path associated with a FileOutput record or a LinIntermediatePath if LinPath points to a workflow and task run subpath.
     * @throws FileNotFoundException if the record does not exist or its type is not a FileOutput or a intermediate directory
     */
    protected Path getTargetOrIntermediatePath() {
        return findTarget(fileSystem, filePath, fragment, false, true)
    }

    /**
     * Get the path associated with a lineage record.
     *
     * @return Path associated with a FileOutput record or a LinMetadataFile with the lineage record for other types, or a intermediate directory
     * @throws FileNotFoundException if the record does not exist
     */
    protected Path getTargetOrMetadataPath() {
        return findTarget(fileSystem, filePath, fragment,true, false)
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
        return "$filePath${query ? '?' + query : ''}${fragment ? '#' + fragment : ''}".toString()
    }
    /**
     * Validates the integrity of the LinPath. If there is a problem with the validation an exception is thrown.
     * To validate just try to get the find target target path. It checks if lid exists, it is a FileOutput,
     * the target path exists and the checksum is the same as the stored in the metadata.
     */
    FileCheck validate() throws Exception{
        final obj = fileSystem.store.load(filePath)
        if( !obj )
            return new FileCheck("File cannot be found")
        if( obj instanceof FileOutput ) {
            final res = validateDataOutput(obj as FileOutput)
            return new FileCheck(res, obj)
        }
        return new FileCheck("Unexpected lineage object type: ${obj.getClass().getName()}")
    }

    @EqualsAndHashCode
    static class FileCheck {
        final String error
        final FileOutput file

        FileCheck(String error, FileOutput out=null) {
            this.error = error
            this.file = out
        }

        /**
         * Implements groovy truth
         */
        boolean asBoolean() {
            return error==null && file!=null
        }
    }
}

