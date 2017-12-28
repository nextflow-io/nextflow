/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
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

package nextflow.file.igfs
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
import org.apache.ignite.igfs.IgfsFile
import org.apache.ignite.igfs.IgfsPath
/**
 * Implements a Ignite file system path compatible with JSR-203
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class IgPath implements Path {

    public static final String PATH_SEPARATOR = '/'

    private IgFileSystem fileSystem

    private Path delegate

    IgPath( IgFileSystem fileSystem, String name ) {
        assert fileSystem
        this.fileSystem = fileSystem
        this.delegate = Paths.get(name)
    }

    IgPath( IgFileSystem fileSystem, String name, String... more ) {
        assert fileSystem
        this.fileSystem = fileSystem
        delegate = Paths.get(name,more)
    }

    private Path unwrap( Path other ) {
        if( other )
            return other instanceof IgPath ? other.delegate : delegate.getFileSystem().getPath(other.toString())
        else
            return null
    }

    @Override
    IgFileSystem getFileSystem() {
        return fileSystem
    }

    @Override
    boolean isAbsolute() {
        return delegate.isAbsolute()
    }

    @Override
    IgPath getRoot() {
        return new IgPath(fileSystem, PATH_SEPARATOR)
    }

    @Override
    IgPath getFileName() {
        return new IgPath(fileSystem, delegate.getFileName().toString())
    }

    @Override
    IgPath getParent() {
        return new IgPath(fileSystem, delegate.getParent().toString())
    }

    @Override
    int getNameCount() {
        return delegate.getNameCount()
    }

    @Override
    IgPath getName(int index) {
        new IgPath(fileSystem, delegate.getName(index).toString())
    }

    @Override
    IgPath subpath(int beginIndex, int endIndex) {
        new IgPath(fileSystem, delegate.subpath(beginIndex, endIndex).toString())
    }

    @Override
    boolean startsWith(Path other) {
        delegate.startsWith(unwrap(other))
    }

    @Override
    boolean startsWith(String other) {
        delegate.startsWith(other)
    }

    @Override
    boolean endsWith(Path other) {
        delegate.endsWith(unwrap(other))
    }

    @Override
    boolean endsWith(String other) {
        delegate.endsWith(other)
    }

    @Override
    IgPath normalize() {
        new IgPath(fileSystem, delegate.normalize().toString())
    }

    @Override
    IgPath resolve(Path other) {
        new IgPath(fileSystem, delegate.resolve(unwrap(other)).toString())
    }

    @Override
    IgPath resolve(String other) {
        new IgPath(fileSystem, delegate.resolve(other).toString())
    }

    @Override
    IgPath resolveSibling(Path other) {
        new IgPath(fileSystem, delegate.resolveSibling(unwrap(other)).toString())
    }

    @Override
    IgPath resolveSibling(String other) {
        new IgPath(fileSystem, delegate.resolveSibling(other).toString())
    }

    @Override
    IgPath relativize(Path other) {
        new IgPath(fileSystem, delegate.relativize(unwrap(other)).toString())
    }

    @Override
    URI toUri() {
        def str = "${IgFileSystemProvider.SCHEME}://${this.delegate.toString()}"
        URI.create(str)
    }

    @Override
    IgPath toAbsolutePath() {
        def str = delegate.toString()
        if( str.startsWith(PATH_SEPARATOR) )
            new IgPath(fileSystem, str)
        else
            new IgPath( fileSystem, '/', str )
    }

    @Override
    IgPath toRealPath(LinkOption... options) throws IOException {
        toAbsolutePath().normalize()
    }

    @Override
    File toFile() {
        throw new UnsupportedOperationException("Method 'toFile' not supported by ${this.class.name}")
    }

    @Memoized
    IgfsPath toIgnitePath() {
        def str = delegate.toString()
        if( str.startsWith(PATH_SEPARATOR) )
            new IgfsPath(str)
        else
            new IgfsPath( '/' + str )
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
    Iterator<IgPath> iterator() {
        return delegate
                .iterator()
                .collect { new IgPath(fileSystem, it.toString()) }
                .iterator()
    }

    @Override
    int compareTo(Path other) {
        return delegate.compareTo(unwrap(other))
    }

    @Override
    String toString() {
        delegate.toString()
    }

    @Override
    boolean equals(Object obj) {
        if (this.is(obj)) return true
        if (getClass() != obj.class) return false

        IgPath other = (IgPath) obj

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
        fileSystem.igfs.delete( toIgnitePath(), recursive )
    }

    void nativeMkdirs() {
        fileSystem.igfs.mkdirs( toIgnitePath() )
    }

    @PackageScope
    IgfsFile nativeReadAttributes() {
        return fileSystem.igfs.info( toIgnitePath() )
    }

    @PackageScope
    def nativeExists() {
        fileSystem.igfs.exists( toIgnitePath() )
    }

    @PackageScope
    def nativeSetTime(long accessTime, long modificationTime) {
        fileSystem.igfs.setTimes(toIgnitePath(), accessTime, modificationTime)
    }
}
