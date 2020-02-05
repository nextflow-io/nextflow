/*
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

package nextflow.processor

import java.nio.file.FileSystem
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.nio.file.Paths
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute

import groovy.transform.CompileStatic
import nextflow.extension.FilesEx
import nextflow.file.FileHolder
import nextflow.util.Escape
import nextflow.util.PathEscapeAware

/**
 * This class model a input file path injected in a process
 * script resolution context and represent the *relative* input file
 * path.
 *
 * This should give access only to metadata information such as name, size, attributes
 * however it's not supposed to allow read/write operations on the file itself
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
final class TaskPath implements Path, PathEscapeAware {

    /**
     * The real path where the file is stored
     */
    private Path target

    /**
     * The relative file name usage to stage it in the
     * task work directory (ie. to symlink it)
     */
    private String alias

    TaskPath(Path source, String alias = null) {
        this.target = source
        this.alias = alias ?: source.getFileName().toString()
    }

    TaskPath(FileHolder holder) {
        this.target = holder.getStorePath()
        this.alias = holder.getStageName()
    }

    // note: required by kryo serialization
    private TaskPath() { }

    @Override
    String toString() {
        return alias
    }

    long size() {
        try {
            return Files.size(target)
        }
        catch( NoSuchFileException e ) {
            return 0
        }
    }

    static boolean equals( Path p1, Path p2 ) {
        if( p1 instanceof TaskPath && p2 instanceof TaskPath ) {
            return p1.equals(p2)
        }
        else if( p1 instanceof TaskPath && p2 ) {
            return p1.alias == p2.getFileName().toString() && p2.equals(p1.target)
        }
        else if( p2 instanceof TaskPath ) {
            return p2.alias == p1.getFileName().toString() && p1.equals(p2.target)
        }
        else {
            return p1?.equals(p2)
        }
    }

    @Override
    boolean equals(Object other) {
        if( other == null )
            return false

        if( other.class != this.class )
            return false

        final that = other as TaskPath
        target.equals(that.target) && alias.equals(that.alias)
    }

    @Override
    int hashCode() {
        Objects.hash(target,alias)
    }

    @Override
    int compareTo(Path other) {
        target.compareTo(other)
    }

    String getName() { alias }

    @Override
    FileSystem getFileSystem() { getFileName().getFileSystem() }

    @Override
    boolean isAbsolute() { false }

    @Override
    Path getRoot() { return null }

    @Override
    Path getFileName() { Paths.get(alias) }

    @Override
    Path getParent() { return null }

    @Override
    int getNameCount() { getFileName().getNameCount() }

    @Override
    Path getName(int index) {
        getFileName().getName(index)
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        getFileName().subpath(beginIndex,endIndex)
    }

    @Override
    boolean startsWith(Path other) {
        return getFileName().startsWith(other.toString())
    }

    @Override
    boolean startsWith(String other) {
        return getFileName().startsWith(other)
    }

    @Override
    boolean endsWith(Path other) {
        return getFileName().endsWith(other.toString())
    }

    @Override
    boolean endsWith(String other) {
        return getFileName().endsWith(other)
    }

    @Override
    Path normalize() {
        return getFileName()
    }

    @Override
    Path resolve(Path other) {
        return getFileName().resolve(other.toString())
    }

    @Override
    Path resolve(String other) {
        return getFileName().resolve(other)
    }

    @Override
    Path resolveSibling(Path other) {
        return getFileName().resolveSibling(other.toString())
    }

    @Override
    Path resolveSibling(String other) {
        return getFileName().resolveSibling(other)
    }

    @Override
    Path relativize(Path other) {
        throw new UnsupportedOperationException()
    }

    @Override
    URI toUri() {
        return getUri()
    }

    @Override
    Path toAbsolutePath() {
        throw new UnsupportedOperationException()
    }

    @Override
    Path toRealPath(LinkOption... options) throws IOException {
        return target
    }

    @Override
    File toFile() {
        getFileName().toFile()
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) throws IOException {
        throw new UnsupportedOperationException()
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>... events) throws IOException {
        throw new UnsupportedOperationException()
    }

    @Override
    Iterator<Path> iterator() {
        return getFileName().iterator()
    }

    String getBaseName() {
        FilesEx.getBaseName(getFileName())
    }

    String getSimpleName() {
        FilesEx.getSimpleName(getFileName())
    }

    String getExtension() {
        FilesEx.getExtension(getFileName())
    }

    boolean exists(LinkOption... options) {
        FilesEx.exists(target,options)
    }

    boolean empty() {
        FilesEx.empty(target)
    }

    boolean canRead() {
        FilesEx.canRead(target)
    }

    boolean canWrite() {
        FilesEx.canWrite(target)
    }

    boolean canExecute() {
        FilesEx.canExecute(target)
    }

    long lastModified(LinkOption...options) {
        FilesEx.lastModified(target,options)
    }

    boolean isHidden() {
        Files.isHidden(getFileName())
    }

    boolean isDirectory(LinkOption... options) {
        Files.isDirectory(target,options)
    }

    boolean isFile(LinkOption...options) {
        Files.isRegularFile(target,options)
    }

    boolean isLink() {
        return true
    }

    Path resolveSymLink() {
        FilesEx.resolveSymLink(target)
    }

    @Deprecated
    BasicFileAttributes readAttributes() {
        FilesEx.readAttributes(target)
    }

    boolean matches( String pattern ) {
        FilesEx.matches(getFileName(), pattern)
    }

    URI getUri() {
        FilesEx.getUri(getFileName())
    }

    Path plus( String other ) {
        FilesEx.plus(getFileName(),other)
    }

    Path plus( Path other ) {
        FilesEx.plus(getFileName(),other)
    }

    Path div( String other ) {
      FilesEx.div(getFileName(),other)
    }

    Path div( Path other ) {
        FilesEx.div(getFileName(),other)
    }

    Path minus(int i) {
        FilesEx.minus(getFileName(),i)
    }

    Path or( String other ) {
        FilesEx.or(getFileName(),other)
    }

    Path or( Path other ) {
        FilesEx.or(getFileName(),other)
    }


    String getPermissions() {
        FilesEx.getPermissions(target)
    }

    boolean setPermissions( String permissions ) {
        throw new UnsupportedOperationException()
    }

    boolean setPermissions( int owner, int group, int other ) {
        throw new UnsupportedOperationException()
    }

    boolean delete() {
        throw new UnsupportedOperationException()
    }

    boolean deleteDir() {
        throw new UnsupportedOperationException()
    }

    void deleteOnExit() {
        throw new UnsupportedOperationException()
    }

    void createDirIfNotExists() {
        throw new UnsupportedOperationException()
    }

    Path copyTo( Path target ) {
        throw new UnsupportedOperationException()
    }

    Path copyTo( String target ) {
        throw new UnsupportedOperationException()
    }

    Path moveTo( Path target ) {
        throw new UnsupportedOperationException()
    }

    Path moveTo( String target ) {
        throw new UnsupportedOperationException()
    }

    boolean mkdir(FileAttribute<?>... attr) {
        throw new UnsupportedOperationException()
    }

    boolean mkdirs(FileAttribute<?>... attr) {
        throw new UnsupportedOperationException()
    }

    boolean renameTo(Path target) {
        throw new UnsupportedOperationException()
    }

    boolean renameTo(String target) {
        throw new UnsupportedOperationException()
    }

    boolean setLastModified(long time) {
        throw new UnsupportedOperationException()
    }

    boolean setExecutable(boolean executable, boolean ownerOnly = true) {
        throw new UnsupportedOperationException()
    }

    boolean setReadable(boolean readable, boolean ownerOnly = true) {
        throw new UnsupportedOperationException()
    }

    boolean setWritable(boolean writable, boolean ownerOnly = true ) {
        throw new UnsupportedOperationException()
    }

    boolean setReadOnly( ) {
        throw new UnsupportedOperationException()
    }

    void rollFile() {
        throw new UnsupportedOperationException()
    }

    Path mklink( Map opts, Path link ) {
        throw new UnsupportedOperationException()
    }

    Path mklink( Map opts, File link ) {
        throw new UnsupportedOperationException()
    }

    Path mklink( Map opts, String link ) {
        throw new UnsupportedOperationException()
    }

    Path mklink( Path link ) {
        throw new UnsupportedOperationException()
    }

    Path mklink( File link ) {
        throw new UnsupportedOperationException()
    }

    Path mklink( String link ) {
        throw new UnsupportedOperationException()
    }

    Path mklink( Map opts = null ) {
        throw new UnsupportedOperationException()
    }

    Path complete() {
        throw new UnsupportedOperationException()
    }

    @Override
    String toStringEscape() {
        return Escape.path(this)
    }
}