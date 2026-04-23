/*
 * Copyright 2013-2026, Seqera Labs
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

package io.seqera.tower.plugin.fs

import java.nio.file.FileSystem
import java.nio.file.InvalidPathException
import java.nio.file.LinkOption
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.WatchEvent
import java.nio.file.WatchKey
import java.nio.file.WatchService

import groovy.transform.CompileStatic

/**
 * {@link Path} implementation for the {@code seqera://} scheme.
 *
 * Path shape:
 * <pre>
 *   seqera://                               depth 0 — root
 *   seqera://&lt;org&gt;                       depth 1
 *   seqera://&lt;org&gt;/&lt;ws&gt;               depth 2
 *   seqera://&lt;org&gt;/&lt;ws&gt;/&lt;type&gt;      depth 3 — resource type
 *   seqera://&lt;org&gt;/&lt;ws&gt;/&lt;type&gt;/...   depth 4+ — handler-owned trail
 * </pre>
 *
 * Resource-type-agnostic for depth &ge; 3: segments after {@code resourceType} are
 * exposed as {@link #getTrail()} and interpreted by the matching
 * {@link ResourceTypeHandler}.
 */
@CompileStatic
class SeqeraPath implements Path {

    public static final String SCHEME = 'seqera'
    public static final String PROTOCOL = "${SCHEME}://"
    public static final String SEPARATOR = '/'

    private final SeqeraFileSystem fs
    private final String org
    private final String workspace
    private final String resourceType
    private final List<String> trail
    /** non-null only for relative paths produced by {@link #relativize(Path)} */
    private final String relPath

    /** Parse a {@code seqera://} URI string. */
    SeqeraPath(SeqeraFileSystem fs, String uriString) {
        this.fs = fs
        this.relPath = null
        if (!uriString.startsWith(PROTOCOL))
            throw new InvalidPathException(uriString, "Not a seqera:// URI")
        final withoutScheme = uriString.substring(PROTOCOL.length())
        final parts = withoutScheme.split('/', -1).toList().findAll { String s -> s != null } as List<String>
        this.org          = parts.size() > 0 && parts[0] ? parts[0] : null
        this.workspace    = parts.size() > 1 && parts[1] ? parts[1] : null
        this.resourceType = parts.size() > 2 && parts[2] ? parts[2] : null
        // Trail: drop any empty segments (trailing slash, accidental double-slashes)
        final List<String> tail = parts.size() > 3
                ? parts.subList(3, parts.size()).findAll { String s -> s } as List<String>
                : new ArrayList<String>()
        this.trail = Collections.unmodifiableList(tail)
        validatePath(uriString)
    }

    /** Programmatic absolute-path constructor. */
    SeqeraPath(SeqeraFileSystem fs, String org, String workspace, String resourceType, List<String> trail) {
        this.fs = fs
        this.relPath = null
        this.org = org
        this.workspace = workspace
        this.resourceType = resourceType
        this.trail = trail != null
                ? Collections.unmodifiableList(new ArrayList<String>(trail))
                : Collections.<String>emptyList()
        validatePath(null)
    }

    /** Relative path produced only by {@link #relativize(Path)}. */
    SeqeraPath(String relPath) {
        this.fs = null
        this.relPath = relPath ?: ''
        this.org = null
        this.workspace = null
        this.resourceType = null
        this.trail = Collections.<String>emptyList()
    }

    private void validatePath(String original) {
        final label = original ?: rawPath()
        if (trail && !resourceType)
            throw new InvalidPathException(label, "Trail segments require a resource-type segment")
        if (resourceType && !workspace)
            throw new InvalidPathException(label, "Resource type requires a workspace segment")
        if (workspace && !org)
            throw new InvalidPathException(label, "Workspace requires an org segment")
        if (org?.contains('/'))
            throw new InvalidPathException(label, "Org name cannot contain '/'")
        if (workspace?.contains('/'))
            throw new InvalidPathException(label, "Workspace name cannot contain '/'")
        if (resourceType?.contains('/'))
            throw new InvalidPathException(label, "Resource type cannot contain '/'")
        for (String t : trail) {
            if (t == null || t.isEmpty())
                throw new InvalidPathException(label, "Path segments cannot be empty")
            if (t.contains('/'))
                throw new InvalidPathException(label, "Path segments cannot contain '/'")
        }
    }

    private String rawPath() {
        final sb = new StringBuilder(PROTOCOL)
        if (org) sb.append(org)
        if (workspace) sb.append('/').append(workspace)
        if (resourceType) sb.append('/').append(resourceType)
        for (String t : trail) sb.append('/').append(t)
        return sb.toString()
    }

    private List<String> nameComponents() {
        if (isAbsolute()) {
            final d = depth()
            final out = new ArrayList<String>(d)
            for (int i = 0; i < d; i++)
                out.add(getName(i).toString())
            return out
        }
        if (!relPath) return Collections.<String>emptyList()
        return relPath.split('/').toList().findAll { String s -> s } as List<String>
    }

    // ---- accessors ----

    String getOrg() { org }
    String getWorkspace() { workspace }
    String getResourceType() { resourceType }
    List<String> getTrail() { trail }

    int depth() {
        if (resourceType) return 3 + trail.size()
        if (workspace) return 2
        if (org) return 1
        return 0
    }

    // ---- Path API ----

    @Override FileSystem getFileSystem() { fs }
    @Override boolean isAbsolute() { fs != null }

    @Override
    Path getRoot() { new SeqeraPath(fs, null, null, null, null) }

    @Override
    Path getFileName() {
        final d = depth()
        if (d == 0) return null
        if (d >= 4) return new SeqeraPath(trail[trail.size() - 1])
        if (d == 3) return new SeqeraPath(resourceType)
        if (d == 2) return new SeqeraPath(workspace)
        return new SeqeraPath(org)
    }

    @Override
    Path getParent() {
        final d = depth()
        if (d == 0) return null
        if (d == 1) return new SeqeraPath(fs, null, null, null, null)
        if (d == 2) return new SeqeraPath(fs, org, null, null, null)
        if (d == 3) return new SeqeraPath(fs, org, workspace, null, null)
        // d >= 4: drop last trail segment
        final newTrail = trail.subList(0, trail.size() - 1)
        return new SeqeraPath(fs, org, workspace, resourceType, newTrail)
    }

    @Override int getNameCount() { depth() }

    @Override
    Path getName(int index) {
        final d = depth()
        if (index < 0 || index >= d)
            throw new IllegalArgumentException("Index out of range: $index")
        if (index == 0) return new SeqeraPath(org)
        if (index == 1) return new SeqeraPath(workspace)
        if (index == 2) return new SeqeraPath(resourceType)
        return new SeqeraPath(trail[index - 3])
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        throw new UnsupportedOperationException("subpath not supported by seqera:// paths")
    }

    @Override
    boolean startsWith(Path other) {
        if (other !instanceof SeqeraPath) return false
        final that = (SeqeraPath) other
        if (this.isAbsolute() != that.isAbsolute()) return false
        final mine = nameComponents()
        final theirs = that.nameComponents()
        if (theirs.size() > mine.size()) return false
        for (int i = 0; i < theirs.size(); i++)
            if (mine[i] != theirs[i]) return false
        return true
    }

    @Override
    boolean startsWith(String other) {
        if (!other) return false
        try {
            final Path p = SeqeraPath.isSeqeraUri(other) ? new SeqeraPath(fs, other) : new SeqeraPath(other)
            return startsWith(p)
        } catch (Exception ignored) { return false }
    }

    @Override
    boolean endsWith(Path other) {
        if (other !instanceof SeqeraPath) return false
        final that = (SeqeraPath) other
        if (that.isAbsolute()) return this.equals(that)
        final mine = nameComponents()
        final theirs = that.nameComponents()
        if (theirs.isEmpty() || theirs.size() > mine.size()) return false
        final offset = mine.size() - theirs.size()
        for (int i = 0; i < theirs.size(); i++)
            if (mine[offset + i] != theirs[i]) return false
        return true
    }

    @Override
    boolean endsWith(String other) {
        if (!other) return false
        try {
            final Path p = SeqeraPath.isSeqeraUri(other) ? new SeqeraPath(fs, other) : new SeqeraPath(other)
            return endsWith(p)
        } catch (Exception ignored) { return false }
    }

    @Override Path normalize() { this }

    @Override
    Path resolve(Path other) {
        if (other instanceof SeqeraPath) {
            final that = (SeqeraPath) other
            if (that.isAbsolute()) return that
            return resolve(that.relPath)
        }
        return resolve(other.toString())
    }

    @Override
    Path resolve(String segment) {
        if (!segment) return this
        if (segment.startsWith(PROTOCOL))
            return new SeqeraPath(fs, segment)
        final stripped = segment.startsWith(SEPARATOR) ? segment.substring(1) : segment
        if (!stripped) return this
        final segs = stripped.split(SEPARATOR, -1).findAll { String s -> s } as List<String>
        SeqeraPath result = this
        for (String seg : segs) result = result.resolveOne(seg)
        return result
    }

    private SeqeraPath resolveOne(String seg) {
        final d = depth()
        if (d == 0) return new SeqeraPath(fs, seg, null, null, null)
        if (d == 1) return new SeqeraPath(fs, org, seg, null, null)
        if (d == 2) return new SeqeraPath(fs, org, workspace, seg, null)
        final newTrail = new ArrayList<String>(trail)
        newTrail.add(seg)
        return new SeqeraPath(fs, org, workspace, resourceType, newTrail)
    }

    @Override
    Path resolveSibling(Path other) {
        final parent = getParent()
        return parent != null ? parent.resolve(other) : other
    }

    @Override
    Path resolveSibling(String other) {
        final parent = getParent()
        return parent != null ? parent.resolve(other) : new SeqeraPath(fs, other)
    }

    @Override
    Path relativize(Path other) {
        if (other !instanceof SeqeraPath) throw new ProviderMismatchException()
        final that = (SeqeraPath) other
        if (!this.isAbsolute() || !that.isAbsolute())
            throw new IllegalArgumentException("Both paths must be absolute to relativize: ${this} vs ${other}")
        final mine = this.nameComponents()
        final theirs = that.nameComponents()
        int common = 0
        while (common < mine.size() && common < theirs.size() && mine[common] == theirs[common]) common++
        final parts = new ArrayList<String>()
        for (int i = common; i < mine.size(); i++) parts.add('..')
        for (int i = common; i < theirs.size(); i++) parts.add(theirs[i])
        return new SeqeraPath(parts.join(SEPARATOR))
    }

    @Override
    URI toUri() {
        String uriPath = null
        if (workspace) {
            final segments = [workspace]
            if (resourceType) segments.add(resourceType)
            for (String t : trail) segments.add(t)
            uriPath = '/' + segments.join('/')
        }
        return new URI(SCHEME, org ?: '', uriPath, null, null)
    }

    @Override
    String toString() {
        if (!isAbsolute()) return relPath
        if (depth() == 0) return PROTOCOL
        return toUri().toString()
    }

    @Override
    Path toAbsolutePath() {
        if (!isAbsolute())
            throw new IllegalStateException("Cannot convert relative SeqeraPath to absolute — no default directory context")
        return this
    }

    @Override Path toRealPath(LinkOption... options) { this }

    @Override
    File toFile() { throw new UnsupportedOperationException("toFile() not supported for seqera:// paths") }

    @Override
    WatchKey register(WatchService w, WatchEvent.Kind<?>[] e, WatchEvent.Modifier... m) {
        throw new UnsupportedOperationException("WatchService not supported by seqera:// paths")
    }

    @Override
    WatchKey register(WatchService w, WatchEvent.Kind<?>... e) {
        throw new UnsupportedOperationException("WatchService not supported by seqera:// paths")
    }

    @Override
    Iterator<Path> iterator() {
        final d = depth()
        final out = new ArrayList<Path>(d)
        for (int i = 0; i < d; i++) out.add(getName(i))
        return out.iterator()
    }

    @Override int compareTo(Path other) { toString().compareTo(other.toString()) }

    @Override
    boolean equals(Object obj) {
        if (obj == this) return true
        if (obj !instanceof SeqeraPath) return false
        return toString() == obj.toString()
    }

    @Override int hashCode() { toString().hashCode() }

    static URI asUri(String path) {
        if (!path) throw new IllegalArgumentException("Missing 'path' argument")
        if (!path.startsWith(PROTOCOL))
            throw new IllegalArgumentException("Invalid Seqera file system path URI - it must start with '${PROTOCOL}' prefix - offending value: $path")
        if (path.startsWith(PROTOCOL + SEPARATOR) && path.length() > PROTOCOL.length() + 1)
            throw new IllegalArgumentException("Invalid Seqera file system path URI - make sure the scheme prefix does not contain more than two slash characters or a query in the root '/' - offending value: $path")
        if (path.startsWith(PROTOCOL + './'))
            path = PROTOCOL + path.substring(PROTOCOL.length() + 2)
        if (path == PROTOCOL || path == PROTOCOL + '.')
            return new URI(PROTOCOL + '/')
        return new URI(path)
    }

    static boolean isSeqeraUri(String path) {
        return path && path.startsWith(PROTOCOL)
    }
}
