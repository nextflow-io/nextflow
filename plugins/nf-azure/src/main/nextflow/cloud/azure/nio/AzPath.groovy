/*
 * Copyright 2021, Microsoft Corp
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
package nextflow.cloud.azure.nio

import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.ProviderMismatchException
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService

import com.azure.storage.blob.BlobClient
import com.azure.storage.blob.BlobContainerClient
import com.azure.storage.blob.models.BlobItem
import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.PackageScope

/**
 * Implements Azure path object
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@EqualsAndHashCode(includes = 'fs,path,directory', includeFields = true)
class AzPath implements Path {

    private AzFileSystem fs

    private Path path

    private AzFileAttributes attributes

    @PackageScope
    boolean directory

    @PackageScope
    AzPath() {}

    @PackageScope
    AzPath( AzFileSystem fs, String path ) {
        this(fs, Paths.get(path), path.endsWith("/") || path=="/$fs.containerName".toString())
    }

    @PackageScope
    AzPath( AzFileSystem fs, BlobClient client ) {
        this(fs, "/${client.getContainerName()}/${client.getBlobName()}")
        this.attributes = new AzFileAttributes(client)
    }

    AzPath(AzFileSystem fs, BlobItem item) {
        this(fs, "/${fs.containerName}/${item.name}")
        this.attributes = new AzFileAttributes(fs.containerName, item)
    }

    @PackageScope
    AzPath(AzFileSystem fs, BlobContainerClient client) {
        this(fs, "/${client.blobContainerName}")
        this.attributes = new AzFileAttributes(client)
    }

    @PackageScope
    AzPath(AzFileSystem fs, Path path, boolean directory) {
        // make sure that the path bucket match the file system bucket
        if( path.isAbsolute() && path.nameCount>0 ) {
            def container = path.getName(0).toString()
            if( container != fs.containerName )
                throw new IllegalArgumentException("Azure path `$container` does not match file system bucket: `${fs.containerName}`")
        }

        this.fs = fs
        this.path = path
        this.directory = directory
    }

    @PackageScope
    AzPath setAttributes(AzFileAttributes attrs) {
        this.attributes = attrs
        return this
    }

    boolean isDirectory() {
        return directory
    }

    String checkContainerName() {
        if( !isAbsolute() )
            throw new IllegalArgumentException("Azure blob container name is not available on relative blob path: $path")
        return path.subpath(0,1)
    }

    BlobContainerClient containerClient() {
        return fs.getBlobServiceClient().getBlobContainerClient(checkContainerName())
    }

    String blobName() {
        if( !path.isAbsolute() )
            return path.toString()

        if( path.nameCount>1 )
            return path.subpath(1, path.nameCount).toString()

        return null
    }

    BlobClient blobClient() {
        def name = blobName()
        if( !name || isContainer() )
            throw new IllegalArgumentException("Azure blob client is not available for container path $path")
        if( directory && !name.endsWith('/') )
            name += '/'
        containerClient().getBlobClient(name)
    }

    @Override
    AzFileSystem getFileSystem() {
        return fs
    }

    @Override
    boolean isAbsolute() {
        path.isAbsolute()
    }

    @Override
    Path getRoot() {
        path.isAbsolute() ? new AzPath(fs, "/${path.getName(0)}/") : null
    }

    @Override
    Path getFileName() {
        final name = path.getFileName()
        name ? new AzPath(fs, name, directory) : null
    }

    @Override
    Path getParent() {
        if( path.isAbsolute() && path.nameCount>1 ) {
            new AzPath(fs, path.parent, true)
        }
        else {
            null
        }
    }

    @Override
    int getNameCount() {
        path.getNameCount()
    }

    @Override
    Path getName(int index) {
        final dir = index < path.getNameCount()-1
        new AzPath(fs, path.getName(index), dir)
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        final dir = endIndex < path.getNameCount()-1
        new AzPath(fs, path.subpath(beginIndex,endIndex), dir)
    }

    @Override
    boolean startsWith(Path other) {
        path.startsWith(other.toString())
    }

    @Override
    boolean startsWith(String other) {
        path.startsWith(other)
    }

    @Override
    boolean endsWith(Path other) {
        path.endsWith(other.toString())
    }

    @Override
    boolean endsWith(String other) {
        path.endsWith(other)
    }

    @Override
    Path normalize() {
        new AzPath(fs, path.normalize(), directory)
    }

    @Override
    AzPath resolve(Path other) {
        if( other.class != AzPath )
            throw new ProviderMismatchException()

        final that = (AzPath)other
        if( other.isAbsolute() )
            return that

        def newPath = path.resolve(that.path)
        new AzPath(fs, newPath, false)
    }

    @Override
    AzPath resolve(String other) {
        if( other.startsWith('/') )
            return (AzPath)fs.provider().getPath(new URI("$AzFileSystemProvider.SCHEME:/$other"))

        def dir = other.endsWith('/')
        def newPath = path.resolve(other)
        new AzPath(fs, newPath, dir)
    }

    @Override
    Path resolveSibling(Path other) {
        if( other.class != AzPath )
            throw new ProviderMismatchException()

        final that = (AzPath)other
        def newPath = path.resolveSibling(that.path)
        if( newPath.isAbsolute() )
            fs.getPath(newPath.toString())
        else
            new AzPath(fs, newPath, false)
    }

    @Override
    Path resolveSibling(String other) {
        def newPath = path.resolveSibling(other)
        if( newPath.isAbsolute() )
            fs.getPath(newPath.toString())
        else
            new AzPath(fs, newPath, false)
    }

    @Override
    Path relativize(Path other) {
        if( other.class != AzPath )
            throw new ProviderMismatchException()

        def newPath = path.relativize( ((AzPath)other).path )
        new AzPath(fs,newPath,false)
    }

    @Override
    String toString() {
        path.toString()
    }

    @Override
    URI toUri() {
        return new URI(toUriString())
    }

    @Override
    Path toAbsolutePath() {
        if(isAbsolute()) return this
        throw new UnsupportedOperationException("Operation 'toAbsolutePath' is not supported by AzPath")
    }

    @Override
    Path toRealPath(LinkOption... options) throws IOException {
        return toAbsolutePath()
    }

    @Override
    File toFile() {
        throw new UnsupportedOperationException("Operation 'toFile' is not supported by AzPath")
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) throws IOException {
        throw new UnsupportedOperationException("Operation 'register' is not supported by AzPath")
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>... events) throws IOException {
        throw new UnsupportedOperationException("Operation 'register' is not supported by AzPath")
    }

    @Override
    Iterator<Path> iterator() {
        final count = path.nameCount
        List<Path> paths = new ArrayList<>()
        for( int i=0; i<count; i++ ) {
            def dir = i<count-1
            paths.add(i, new AzPath(fs, path.getName(i), dir))
        }
        paths.iterator()
    }

    @Override
    int compareTo(Path other) {
        return this.toString() <=> other.toString()
    }

    String getContainerName() {
        if( path.isAbsolute() ) {
            path.nameCount==0 ? '/' : path.getName(0)
        }
        else
            return null
    }

    boolean isContainer() {
        path.isAbsolute() && path.nameCount<2
    }


    String toUriString() {
        if( path.isAbsolute() ) {
            return "${AzFileSystemProvider.SCHEME}:/${path.toString()}"
        }
        else {
            return "${AzFileSystemProvider.SCHEME}:${path.toString()}"
        }
    }

    AzFileAttributes attributesCache() {
        def result = attributes
        attributes = null
        return result
    }

}
