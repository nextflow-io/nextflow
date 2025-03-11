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
 *
 */

package nextflow.data.cid.fs

import groovy.json.JsonSlurper
import groovy.util.logging.Slf4j
import nextflow.data.cid.model.DataType
import nextflow.util.CacheHelper
import nextflow.util.TestOnly

import static nextflow.data.cid.fs.CidFileSystemProvider.*

import java.nio.file.FileSystem
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService

import groovy.transform.CompileStatic
import nextflow.file.FileHelper

/**
 * CID file system path
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class CidPath implements Path {

    static public final String SEPARATOR = '/'
    public static final String CID_PROT = "${SCHEME}://"

    static private final String[] EMPTY = new String[] {}

    private CidFileSystem fileSystem

    // String with the cid file path
    private String filePath

    /*
     * Only needed to prevent serialization issues - see https://github.com/nextflow-io/nextflow/issues/5208
     */
    protected CidPath(){}

    CidPath(CidFileSystem fs, String path) {
        this(fs, path, EMPTY)
    }

    CidPath(CidFileSystem fs, String path, String[] more) {
        this.fileSystem = fs
        this.filePath = resolve0(fs, norm0(path), norm0(more))
    }

    private static void validateHash(Map cidObject) {
        final hashedPath = FileHelper.toCanonicalPath(cidObject.path as String)
        if( !hashedPath.exists() )
            throw new FileNotFoundException("Target path $cidObject.path does not exists.")
        if( cidObject.checksum && CacheHelper.hasher(hashedPath).hash().toString() != cidObject.checksum ) {
            log.warn("Checksum of $hashedPath does not match with the one stored in the metadata")
        }
    }

    @TestOnly
    protected String getFilePath(){ this.filePath }


    /**
     * Finds the target path of a CID path
     **/
    protected static Path findTarget(CidFileSystem fs, String filePath, String[] childs=[]) throws Exception{
        if( !fs )
            throw new IllegalArgumentException("Cannot get target path for a relative CidPath")
        if( filePath.isEmpty() || filePath == SEPARATOR )
            throw new IllegalArgumentException("Cannot get target path for an empty CidPath")
        final store = fs.getCidStore()
        if( !store )
            throw new Exception("CID store not found. Check Nextflow configuration.")
        final slurper = new JsonSlurper()
        final object = store.load(filePath)
        if ( object ){
            final cidObject = slurper.parse(object.toString().toCharArray()) as Map
            final type = DataType.valueOf(cidObject.type as String)
            if( type == DataType.TaskOutput || type == DataType.WorkflowOutput ) {
                // return the real path stored in the metadata
                validateHash(cidObject)
                def realPath = FileHelper.toCanonicalPath(cidObject.path as String)
                if (childs && childs.size() > 0)
                    realPath = realPath.resolve(childs.join(SEPARATOR))
                if( !realPath.exists() )
                    throw new FileNotFoundException("Target path $realPath for $filePath does not exists.")
                return realPath
            }
        } else {
            // If there isn't metadata check the parent to check if it is a subfolder of a task/workflow output
            final currentPath = Path.of(filePath)
            final parent = Path.of(filePath).getParent()
            if( parent) {
                ArrayList<String> newChilds = new ArrayList<String>()
                newChilds.add(currentPath.getFileName().toString())
                newChilds.addAll(childs)
                return findTarget(fs, parent.toString(), newChilds as String[])
            }
        }
        throw new FileNotFoundException("Target path $filePath does not exists.")
    }

    private static boolean isEmptyBase(CidFileSystem fs, String base){
        return !base || base == SEPARATOR || (fs && base == "..")
    }

    private static String resolve0(CidFileSystem fs, String base, String[] more) {
        if( isEmptyBase(fs,base) ) {
            return resolveEmptyPathCase(fs, more as List)
        }
        if( base.contains(SEPARATOR) ) {
            final parts = base.tokenize(SEPARATOR)
            final remain = parts[1..-1] + more.toList()
            return resolve0(fs, parts[0], remain as String[])
        }
        def result = Path.of(base)
        return more ? result.resolve(more.join(SEPARATOR)).toString() : result.toString()
    }

    private static String resolveEmptyPathCase(CidFileSystem fs, List<String> more ){
        switch(more.size()) {
            case 0:
                return "/"
            case 1:
                return resolve0(fs, more[0], EMPTY)
            default:
                return resolve0(fs, more[0], more[1..-1] as String[])
        }
    }

    static private String norm0(String path) {
        if( !path || path==SEPARATOR)
            return ""
        //Remove repeated elements
        path = Path.of(path).normalize().toString()
        //Remove initial and final separators
        if( path.startsWith(SEPARATOR) )
            path = path.substring(1)
        if( path.endsWith(SEPARATOR) )
            path = path.substring(0,path.size()-1)
        return path
    }
    
    static private String[] norm0(String... path) {
        for( int i=0; i<path.length; i++ ) {
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
        return new CidPath(fileSystem, SEPARATOR)
    }

    @Override
    Path getFileName() {
        final result = Path.of(filePath).getFileName()?.toString()
        return result ? new CidPath(null, result) : null
    }

    @Override
    Path getParent() {
        final c = getNameCount()
        if( c>1 )
            return subpath(0,c-1)
        if( c==1 )
            return new CidPath(fileSystem,SEPARATOR)
        return null
    }

    @Override
    int getNameCount() {
        return Path.of(filePath).nameCount
    }

    @Override
    Path getName(int index) {
        if( index<0 )
            throw new IllegalArgumentException("Path name index cannot be less than zero - offending value: $index")
        final path = Path.of(filePath)
        return new CidPath(index==0 ? fileSystem : null, path.getName(index).toString())
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        if( beginIndex<0 )
            throw new IllegalArgumentException("subpath begin index cannot be less than zero - offending value: $beginIndex")
        final path = Path.of(filePath)
        return new CidPath(beginIndex==0 ? fileSystem : null, path.subpath(beginIndex, endIndex).toString())
    }

    @Override
    Path normalize() {
        return new CidPath(fileSystem, Path.of(filePath).normalize().toString())
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
        if( CidPath.class != other.class )
            throw new ProviderMismatchException()

        final that = (CidPath)other

        if( that.fileSystem && this.fileSystem != that.fileSystem )
            return other
        if( that.isAbsolute() ) {
            return that
        } else {
            final newPath = Path.of(filePath).resolve(that.toString())
            return new CidPath(fileSystem, newPath.toString())
        }
    }

    @Override
    Path resolve(String path) {
        if( !path )
            return this
        final scheme = FileHelper.getUrlProtocol(path)
        if( !scheme ) {
            // consider the path as a cid relative path
            return resolve(new CidPath(null,path))
        }
        if( scheme != SCHEME ) {
            throw new ProviderMismatchException()
        }
        final that = fileSystem.provider().getPath(asUri(path))
        return resolve(that)
    }

    @Override
    Path relativize(Path other) {
        if( CidPath.class != other.class ) {
            throw new ProviderMismatchException()
        }
        CidPath cidOther = other as CidPath
        if( this.isAbsolute() != cidOther.isAbsolute() )
            throw new IllegalArgumentException("Cannot compare absolute with relative paths");
        def path
        if( this.isAbsolute() ) {
            // Compare 'filePath' as absolute paths adding the root separator
            path = Path.of(SEPARATOR + filePath).relativize(Path.of(SEPARATOR + cidOther.filePath))
        } else {
            // Compare 'filePath' as relative paths
            path = Path.of(filePath).relativize(Path.of(cidOther.filePath))
        }
        return new CidPath(null , path.getNameCount()>0 ? path.toString(): SEPARATOR)
    }

    @Override
    URI toUri() {
        asUri("${SCHEME}://${filePath}")
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
        return this
    }

    protected Path getTargetPath(){
        return findTarget(fileSystem, filePath)
    }

    @Override
    File toFile() throws IOException {
        throw new UnsupportedOperationException("toFile not supported by CidPath")
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) throws IOException {
        throw new UnsupportedOperationException("Register not supported by CidPath")
    }

    @Override
    int compareTo(Path other) {
        if( CidPath.class != other.class )
            throw new ProviderMismatchException()
        final that = other as CidPath
        return Path.of(this.filePath).compareTo(Path.of(that.filePath))
    }

    @Override
    boolean equals(Object other) {
        if( CidPath.class != other.class ) {
            return false
        }
        final that = (CidPath)other
        return this.fileSystem == that.fileSystem && this.filePath.equals(that.filePath)
    }

    /**
     * @return The unique hash code for this path
     */
    @Override
    int hashCode() {
        return Objects.hash(fileSystem,filePath)
    }

    static URI asUri(String path) {
        if (!path)
            throw new IllegalArgumentException("Missing 'path' argument")
        if (!path.startsWith(CID_PROT))
            throw new IllegalArgumentException("Invalid CID file system path URI - it must start with '${CID_PROT}' prefix - offendinf value: $path")
        if (path.startsWith(CID_PROT + SEPARATOR) && path.length() > 7)
            throw new IllegalArgumentException("Invalid CID file system path URI - make sure the schema prefix does not container more than two slash characters - offending value: $path")
        if (path == CID_PROT) //Empty path case
            return new URI("")
        return new URI(path)
    }

    @Override
    String toString() {
        filePath
    }

}
