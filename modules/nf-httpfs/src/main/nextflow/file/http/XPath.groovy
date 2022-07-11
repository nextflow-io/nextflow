/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

package nextflow.file.http

import java.nio.file.FileSystem
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.ProviderMismatchException
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService

import groovy.transform.CompileStatic
/**
 * Implements a {@link Path} for http/ftp protocols
 *
 * @author Emilio Palumbo <emilio.palumbo@crg.eu>
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class XPath implements Path {

    private static final String[] EMPTY = []

    private XFileSystem fs

    private Path path

    private String query

    XPath(XFileSystem fs, String path) {
        this(fs, path, EMPTY)
    }

    XPath(XFileSystem fs, String path, String[] more) {
        this.fs = fs
        this.query = query(path)
        this.path = Paths.get(stripQuery(path) ?:'/', more)
    }

    private XPath(XFileSystem fs, Path path, String query=null) {
        this.fs = fs
        this.path = path
        this.query = query
    }

    private URI getBaseUri() {
        fs?.getBaseUri()
    }

    private XPath createXPath(String path) {
        fs && path.startsWith('/') ? new XPath(fs, path) : new XPath(null, path)
    }

    @Override
    FileSystem getFileSystem() {
        return fs
    }

    @Override
    boolean isAbsolute() {
        return path.isAbsolute()
    }

    @Override
    Path getRoot() {
        return createXPath("/")
    }

    @Override
    Path getFileName() {
        final result = path?.getFileName()?.toString()
        return result ? new XPath(null, result) : null
    }

    @Override
    Path getParent() {
        String result = path.parent ? path.parent.toString() : null
        if( result ) {
            if( result != '/' ) result += '/'
            return createXPath(result)
        }
        return null
    }

    @Override
    int getNameCount() {
        return path.toString() ? path.nameCount : 0
    }

    @Override
    Path getName(int index) {
        return new XPath(null, path.getName(index).toString())
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        return new XPath(null, path.subpath(beginIndex, endIndex).toString())
    }

    @Override
    Path normalize() {
        return new XPath(fs, path.normalize(), query)
    }

    @Override
    boolean startsWith(Path other) {
        return startsWith(other.toString())
    }

    @Override
    boolean startsWith(String other) {
        return path.startsWith(other)
    }

    @Override
    boolean endsWith(Path other) {
        return endsWith(other.toString())
    }

    @Override
    boolean endsWith(String other) {
        return path.endsWith(other)
    }

    @Override
    Path resolve(Path other) {
        if( this.class != other.class )
            throw new ProviderMismatchException()

        def that = (XPath)other

        if( that.fs && this.fs != that.fs )
            return other

        else if( that.path ) {
            def newPath = this.path.resolve(that.path)
            return new XPath(fs, newPath, that.query)
        }
        else {
            return this
        }

    }

    @Override
    Path resolve(String other) {
        resolve(get(other))
    }

    @Override
    Path resolveSibling(Path other) {
        if( this.class != other.class )
            throw new ProviderMismatchException()

        def that = (XPath)other

        if( that.fs && this.fs != that.fs )
            return other

        if( that.path ) {
            final Path newPath = this.path.resolveSibling(that.path)
            return newPath.isAbsolute() ? new XPath(fs, newPath) : new XPath(null, newPath)
        }
        else {
            return this
        }
    }

    @Override
    Path resolveSibling(String other) {
        resolveSibling(get(other))
    }

    @Override
    Path relativize(Path other) {
        def otherPath = ((XPath)other).path
        return createXPath(path.relativize(otherPath).toString())
    }

    @Override
    URI toUri() {
        final String concat = query!=null ? "$path?$query" : path.toString()
        return baseUri ? new URI("$baseUri$concat") : new URI(concat)
    }

    @Override
    Path toAbsolutePath() {
        return this
    }

    @Override
    Path toRealPath(LinkOption... options) throws IOException {
        return this
    }

    @Override
    File toFile() {
        throw new UnsupportedOperationException()
    }

    @Override
    String toString() {
        return path.toString()
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) throws IOException {
        throw new UnsupportedOperationException("Register not supported by XFileSystem")
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>... events) throws IOException {
        throw new UnsupportedOperationException("Register not supported by XFileSystem")
    }

    @Override
    Iterator<Path> iterator() {
        final len = getNameCount()
        new Iterator<Path>() {
            int index
            Path current = len ? getName(index++) : null

            @Override
            boolean hasNext() {
                return current != null
            }

            @Override
            Path next() {
                final result = current
                current = index<len ? getName(index++) : null
                return result
            }

            @Override
            void remove() {
                throw new UnsupportedOperationException("Remove operation not supported")
            }
        }
    }

    @Override
    int compareTo(Path other) {
        return this.toUri().toString() <=> other.toUri().toString()
    }

    @Override
    boolean equals(Object other) {
        if (other.class != XPath) {
            return false
        }
        final that = (XPath)other
        return this.fs == that.fs && this.path == that.path && this.query == that.query
    }

    /**
     * @return The unique hash code for this pah
     */
    @Override
    int hashCode() {
        return Objects.hash(fs,path,query)
    }

    /**
     * Path factory method.
     *
     * @param str
     *      A fully qualified URI path string e.g. {@code http://www.host.name/data.file.txt}
     * @return
     *      The corresponding {@link XPath} object
     * @throws
     *      ProviderMismatchException When the URI scheme does not match any supported protocol i.e. {@code ftp},
     *      {@code http}, {@code https}
     */
    static XPath get(String str) {
        if( str == null )
            return null

        def uri = new URI(null,null,str,null,null)

        if( uri.scheme && !XFileSystemProvider.ALL_SCHEMES.contains(uri.scheme))
            throw new ProviderMismatchException()

        uri.authority ? (XPath)Paths.get(uri) : new XPath(null, str)
    }

    static String stripQuery(String uri) {
        if(!uri)
            return null
        final p = uri.indexOf('?')
        return p>0 ? uri.substring(0,p) : (p==0 ? null : uri)
    }

    static String query(String uri) {
        if(!uri)
            return null
        final p = uri.indexOf('?')
        return p!=-1 && p<uri.size()-1 ? uri.substring(p+1) : null
    }
}
