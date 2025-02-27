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

    static public String SEPARATOR = '/'
    public static final String METADATA_FILE = '.data.json'
    public static final String CID_PROT = "${SCHEME}://".toString()

    static private String[] EMPTY = new String[] {}

    private CidFileSystem fileSystem

    // Path of the file in the metadata cid store
    private Path storePath

    // String with the cid file path
    private String filePath

    /*
     * Only needed to prevent serialization issues - see https://github.com/nextflow-io/nextflow/issues/5208
     */
    protected CidPath(){}

    protected CidPath(CidFileSystem fs, Path target) {
        this.fileSystem = fs
        this.storePath = target
        this.filePath = filePath0(fs, target)
    }

    CidPath(CidFileSystem fs, String path) {
        this(fs, path, EMPTY)
    }

    CidPath(CidFileSystem fs, String path, String[] more) {
        this.fileSystem = fs
        this.storePath = resolve0(fs, norm0(path), norm0(more))
        this.filePath = filePath0(fs, storePath)
    }

    @TestOnly
    protected String getFilePath(){ this.filePath }

    @TestOnly
    protected Path getStorePath(){ this.storePath }


    /**
     * Finds the target path of a CID path
     **/
    protected static Path findTarget(Path cidStorePath, CidFileSystem fs, String[] childs=[]){
        assert fs
        if( fs.basePath == cidStorePath )
            return null
        final metadata = cidStorePath.resolve(METADATA_FILE).toFile()
        if ( metadata.exists() ){
            final slurper = new JsonSlurper()
            final cidObject = slurper.parse(metadata.text.toCharArray()) as Map
            final type = DataType.valueOf(cidObject.type as String)
            if( type == DataType.TaskOutput || type == DataType.WorkflowOutput ) {
                // return the real path stored in the metadata
                final realPath = Path.of(cidObject.path as String, childs)
                if( !realPath.exists() )
                    throw new FileNotFoundException("Target path $realPath for $cidStorePath does not exists.")
                if( cidObject.checksum && CacheHelper.hasher(realPath).hash().toString() != cidObject.checksum ) {
                    log.warn("Checksum of $cidStorePath does not match with the one stored in the metadata")
                }
                return realPath
            }
        } else {
            // If there isn't metadata check the parent to check if it is a subfolder of a task/workflow output
            final parent = cidStorePath.getParent()
            if( parent) {
                ArrayList<String> newChilds = new ArrayList<String>()
                newChilds.add(cidStorePath.getFileName().toString())
                newChilds.addAll(childs)
                return findTarget(parent, fs, newChilds as String[])
            }
        }
        return null
    }

    private static String filePath0(CidFileSystem fs, Path target) {
        if( !fs )
            return target.toString()
        return fs.basePath != target
                ? fs.basePath.relativize(target).toString()
                : SEPARATOR
    }

    private static Path resolve0(CidFileSystem fs, String base, String[] more) {
        if( !base || base == SEPARATOR ) {
            return resolveEmptyPathCase(fs, more as List)
        }
        if( base.contains(SEPARATOR) ) {
            final parts = base.tokenize(SEPARATOR)
            final remain = parts[1..-1] + more.toList()
            return resolve0(fs, parts[0], remain as String[])
        }
        final result = fs ? fs.basePath.resolve(base) : Path.of(base)
        return more
            ? result.resolve(more.join(SEPARATOR))
            : result
    }

    private static Path resolveEmptyPathCase(CidFileSystem fs, List<String> more ){
        switch(more.size()) {
            case 0:
                return fs ? fs.basePath : Path.of("/")
            case 1:
                return resolve0(fs, more[0], EMPTY)
            default:
                return resolve0(fs, more[0], more[1..-1] as String[])
        }

    }

    static private String norm0(String path) {
        if( !path )
            return ""
        if( path==SEPARATOR )
            return path
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
        final result = storePath?.getFileName()?.toString()
        return result ? new CidPath(null, result) : null
    }

    @Override
    Path getParent() {
        final c = getNameCount()
        if( c>1 )
            return subpath(0,c-1)
        if( c==1 )
            return new CidPath(fileSystem,"/")
        return null
    }

    @Override
    int getNameCount() {
        return fileSystem ? storePath.nameCount-fileSystem.basePath.nameCount : storePath.nameCount
    }

    @Override
    Path getName(int index) {
        if( index<0 )
            throw new IllegalArgumentException("Path name index cannot be less than zero - offending value: $index")
        final c= fileSystem.basePath.nameCount
        return new CidPath(index==0 ? fileSystem : null, storePath.getName(c + index).toString())
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        if( beginIndex<0 )
            throw new IllegalArgumentException("subpath begin index cannot be less than zero - offending value: $beginIndex")
        final c= fileSystem.basePath.nameCount
        return new CidPath(beginIndex==0 ? fileSystem : null, storePath.subpath(c+beginIndex, c+endIndex).toString())
    }

    @Override
    Path normalize() {
        return new CidPath(fileSystem, storePath.normalize())
    }

    @Override
    boolean startsWith(Path other) {
        return startsWith(other.toString())
    }

    @Override
    boolean startsWith(String other) {
        return storePath.startsWith(fileSystem.basePath.resolve(other))
    }

    @Override
    boolean endsWith(Path other) {
        return endsWith(other.toString())
    }

    @Override
    boolean endsWith(String other) {
        return storePath.endsWith(other)
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
        }
        if( that.storePath ) {
            final newPath = this.storePath.resolve(that.storePath)
            return new CidPath(fileSystem, newPath)
        }
        return this
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
        final path = storePath.relativize(((CidPath) other).storePath)
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
        return getTargetPath()
    }

    protected Path getTargetPath(){
        final target = findTarget(storePath, fileSystem)
        return target ? target : storePath
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
        return this.storePath.compareTo(that.storePath)
    }

    @Override
    boolean equals(Object other) {
        if( CidPath.class != other.class ) {
            return false
        }
        final that = (CidPath)other
        return this.fileSystem == that.fileSystem && this.storePath.equals(that.storePath)
    }

    /**
     * @return The unique hash code for this path
     */
    @Override
    int hashCode() {
        return Objects.hash(fileSystem,storePath)
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
