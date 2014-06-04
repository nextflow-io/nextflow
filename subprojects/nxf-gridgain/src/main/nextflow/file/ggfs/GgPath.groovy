/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.file.ggfs
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import org.gridgain.grid.ggfs.GridGgfsFile
import org.gridgain.grid.ggfs.GridGgfsPath

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class GgPath implements Path {

    public static final String PATH_SEPARATOR = '/'

    private GgFileSystem fileSystem

    private Path delegate

    GgPath( GgFileSystem fileSystem, String name ) {
        assert fileSystem
        this.fileSystem = fileSystem
        this.delegate = Paths.get(name)
    }

    GgPath( GgFileSystem fileSystem, String name, String... more ) {
        assert fileSystem
        this.fileSystem = fileSystem
        delegate = Paths.get(name,more)
    }

    @Override
    GgFileSystem getFileSystem() {
        return fileSystem
    }

    @Override
    boolean isAbsolute() {
        return delegate.isAbsolute()
    }

    @Override
    GgPath getRoot() {
        return new GgPath(fileSystem, PATH_SEPARATOR)
    }

    @Override
    GgPath getFileName() {
        return new GgPath(fileSystem, delegate.getFileName().toString())
    }

    @Override
    GgPath getParent() {
        return new GgPath(fileSystem, delegate.getParent().toString())
    }

    @Override
    int getNameCount() {
        return delegate.getNameCount()
    }

    @Override
    GgPath getName(int index) {
        new GgPath(fileSystem, delegate.getName(index).toString())
    }

    @Override
    GgPath subpath(int beginIndex, int endIndex) {
        new GgPath(fileSystem, delegate.subpath(beginIndex, endIndex).toString())
    }

    @Override
    boolean startsWith(Path other) {
        delegate.startsWith((other as GgPath).delegate)
    }

    @Override
    boolean startsWith(String other) {
        delegate.startsWith(other)
    }

    @Override
    boolean endsWith(Path other) {
        delegate.endsWith((other as GgPath).delegate)
    }

    @Override
    boolean endsWith(String other) {
        delegate.endsWith(other)
    }

    @Override
    GgPath normalize() {
        new GgPath(fileSystem, delegate.normalize().toString())
    }

    @Override
    GgPath resolve(Path other) {
        new GgPath(fileSystem, delegate.resolve((other as GgPath).delegate).toString())
    }

    @Override
    GgPath resolve(String other) {
        new GgPath(fileSystem, delegate.resolve(other).toString())
    }

    @Override
    GgPath resolveSibling(Path other) {
        new GgPath(fileSystem, delegate.resolveSibling((other as GgPath).delegate).toString())
    }

    @Override
    GgPath resolveSibling(String other) {
        new GgPath(fileSystem, delegate.resolveSibling(other).toString())
    }

    @Override
    GgPath relativize(Path other) {
        new GgPath(fileSystem, delegate.relativize((other as GgPath).delegate).toString())
    }

    @Override
    URI toUri() {
        def str = "${GgFileSystemProvider.SCHEME}://${this.delegate.toString()}"
        URI.create(str)
    }

    @Override
    GgPath toAbsolutePath() {
        def str = delegate.toString()
        if( str.startsWith(PATH_SEPARATOR) )
            new GgPath(fileSystem, str)
        else
            new GgPath( fileSystem, '/', str )
    }

    @Override
    GgPath toRealPath(LinkOption... options) throws IOException {
        toAbsolutePath().normalize()
    }

    @Override
    File toFile() {
        throw new UnsupportedOperationException("Method 'toFile' not supported by ${this.class.name}")
    }

    @Memoized
    GridGgfsPath toGridGgfsPath() {
        def str = delegate.toString()
        if( str.startsWith(PATH_SEPARATOR) )
            new GridGgfsPath(str)
        else
            new GridGgfsPath( '/' + str )
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) throws IOException {
        throw new UnsupportedOperationException("Method 'register' not supported by ${this.class.name}")
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>... events) throws IOException {
        throw new UnsupportedOperationException("Method 'register' not supported by ${this.class.name}")
    }

    @Override
    Iterator<GgPath> iterator() {
        return delegate
                .iterator()
                .collect { new GgPath(fileSystem, it.toString()) }
                .iterator()
    }

    @Override
    int compareTo(Path other) {
        return delegate.compareTo( (other as GgPath).delegate )
    }

    @Override
    String toString() {
        delegate.toString()
    }

    @Override
    boolean equals(Object obj) {
        if (this.is(obj)) return true
        if (getClass() != obj.class) return false

        GgPath other = (GgPath) obj

        if (!fileSystem.is(other.fileSystem) ) return false
        if (delegate != other.delegate) return false

        return true
    }

    int hashCode() {
        int result
        result = fileSystem.hashCode()
        result = 31 * result + delegate.hashCode()
        return result
    }

    @PackageScope()
    boolean nativeDelete( boolean recursive = false ) {
        fileSystem.ggfs.delete( toGridGgfsPath(), recursive )
    }

    void nativeMkdirs() {
        fileSystem.ggfs.mkdirs( toGridGgfsPath() )
    }

    @PackageScope
    GridGgfsFile nativeReadAttributes() {
        return fileSystem.ggfs.info( toGridGgfsPath() )
    }

    @PackageScope
    def nativeExists() {
        fileSystem.ggfs.exists( toGridGgfsPath() )
    }

    @PackageScope
    def nativeSetTime(long accessTime, long modificationTime) {
        fileSystem.ggfs.setTimes(toGridGgfsPath(), accessTime, modificationTime)
    }
}
