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

package nextflow.data.fs

import static nextflow.data.fs.CidFileSystemProvider.*

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
 * Model a CID file system path
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CidPath implements Path {

    static public String SEPARATOR = '/'

    static private String[] EMTPY = new String[] {}

    private CidFileSystem fileSystem

    private Path target

    private String filePath

    /*
     * Only needed to prevent serialization issues - see https://github.com/nextflow-io/nextflow/issues/5208
     */
    protected CidPath(){}

    protected CidPath(CidFileSystem fs, Path target) {
        this.fileSystem = fs
        this.target = target
        this.filePath = filePath0(fs, target)
    }

    CidPath(CidFileSystem fs, String path) {
        this(fs, norm0(path), EMTPY)
    }

    CidPath(CidFileSystem fs, String path, String[] more) {
        this.fileSystem = fs
        this.target = resolve0(fs, norm0(path), norm0(more))
        this.filePath = filePath0(fs, target)
    }

    private static String filePath0(CidFileSystem fs, Path target) {
        if( !fs )
            return target.toString()
        return fs.basePath != target
                ? fs.basePath.relativize(target).toString()
                : SEPARATOR
    }

    private static Path resolve0(CidFileSystem fs, String base, String[] more) {
        if( !base )
            throw new IllegalArgumentException("Missing CID base path")
        if( base==SEPARATOR && !more )
            return fs.basePath
        if( base.contains(SEPARATOR) ) {
            final parts = base.tokenize(SEPARATOR)
            final remain = parts[1..-1] + more.toList()
            resolve0(fs, parts[0], remain as String[])
        }
        final result = fs ? fs.basePath.resolve(base) : Path.of(base)
        return more
            ? result.resolve(more.join(SEPARATOR))
            : result
    }

    static private String norm0(String path) {
        if( path==SEPARATOR )
            return path
        if( !path )
            return path
        while( path.startsWith('/') )
            path = path.substring(1)
        while( path.endsWith('/') )
            path = path.substring(0,path.size()-1)
        return path
    }
    
    static private String[] norm0(String... path) {
        for( int i=0; i<path.length; i++ ) {
            path[i] = norm0(path[i])
        }
        return path
    }

    CidFileId getFileId() {
        return new CidFileId(filePath)
    }

    @Override
    FileSystem getFileSystem() {
        return fileSystem
    }

    FileSystem getTargetSystem() {
        return target.getFileSystem()
    }

    @Override
    boolean isAbsolute() {
        return fileSystem!=null
    }

    @Override
    Path getRoot() {
        return new CidPath(fileSystem, SEPARATOR)
    }

    @Override
    Path getFileName() {
        final result = target?.getFileName()?.toString()
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
        return target.toString() ? target.nameCount-fileSystem.basePath.nameCount : 0
    }

    @Override
    Path getName(int index) {
        if( index<0 )
            throw new IllegalArgumentException("Path name index cannot be less than zero - offending value: $index")
        final c= fileSystem.basePath.nameCount
        return new CidPath(index==0 ? fileSystem : null, target.getName(c + index).toString())
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        if( beginIndex<0 )
            throw new IllegalArgumentException("subpath begin index cannot be less than zero - offending value: $beginIndex")
        final c= fileSystem.basePath.nameCount
        return new CidPath(beginIndex==0 ? fileSystem : null, target.subpath(c+beginIndex, c+endIndex).toString())
    }

    @Override
    Path normalize() {
        return new CidPath(fileSystem, target.normalize())
    }

    @Override
    boolean startsWith(Path other) {
        return startsWith(other.toString())
    }

    @Override
    boolean startsWith(String other) {
        return target.startsWith(fileSystem.basePath.resolve(other))
    }

    @Override
    boolean endsWith(Path other) {
        return endsWith(other.toString())
    }

    @Override
    boolean endsWith(String other) {
        return target.endsWith(other)
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
        if( that.target ) {
            final newPath = this.target.resolve(that.target)
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
        final otherPath = ((CidPath)other).target
        return new CidPath(fileSystem, target.relativize(otherPath).toString())
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
        return target
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
        return this.target.compareTo(that.target)
    }

    @Override
    boolean equals(Object other) {
        if( CidPath.class != other.class ) {
            return false
        }
        final that = (CidPath)other
        return this.fileSystem == that.fileSystem && this.target == that.target
    }

    /**
     * @return The unique hash code for this path
     */
    @Override
    int hashCode() {
        return Objects.hash(fileSystem,target)
    }

    static URI asUri(String path) {
        if( !path )
            throw new IllegalArgumentException("Missing 'path' argument")
        if( !path.startsWith('cid://') )
            throw new IllegalArgumentException("Invalid CID file system path URI - it must start with 'cid://' prefix - offendinf value: $path")
        if( path.startsWith('cid:///') && path.length()>7 )
            throw new IllegalArgumentException("Invalid CID file system path URI - make sure the schema prefix does not container more than two slash characters - offending value: $path")

        // note: this URI constructor parse the path parameter and extract the `scheme` and `authority` components
        return new URI(null,null, path,null,null)
    }

    @Override
    String toString() {
        return toUriString()
    }
}
