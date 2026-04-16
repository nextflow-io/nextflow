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
 * Path hierarchy:
 * <pre>
 *   depth 0  seqera://                                   (root — directory)
 *   depth 1  seqera://&lt;org&gt;                              (org — directory)
 *   depth 2  seqera://&lt;org&gt;/&lt;workspace&gt;                  (workspace — directory)
 *   depth 3  seqera://&lt;org&gt;/&lt;workspace&gt;/datasets          (resource type — directory)
 *   depth 4  seqera://&lt;org&gt;/&lt;workspace&gt;/datasets/&lt;name&gt;  (dataset file)
 *            seqera://&lt;org&gt;/&lt;workspace&gt;/datasets/&lt;name@ver&gt;  (pinned version)
 * </pre>
 *
 * @author Seqera Labs
 */
@CompileStatic
class SeqeraPath implements Path {

    /** URI scheme */
    public static final String SCHEME = 'seqera'
    public static final String PROTOCOL = "${SCHEME}://"
    public static final String SEPARATOR = '/'

    private final SeqeraFileSystem fs
    /** path segments in order: [org, workspace, resourceType, datasetName] — null for missing levels */
    private final String org
    private final String workspace
    private final String resourceType
    private final String datasetName
    /** version string extracted from {@code @version} suffix; null when not pinned */
    private final String version
    /**
     * Raw relative path string — non-null only for relative {@code SeqeraPath} instances
     * created by {@link #relativize(Path)}. When non-null, {@link #fs} is {@code null}
     * and all segment fields are {@code null}.
     */
    private final String relPath

    /**
     * Parse a {@code seqera://} URI string into a SeqeraPath.
     * The URI authority is the org; path segments are workspace, resourceType, datasetName.
     * The last segment may contain a {@code @version} suffix.
     */
    SeqeraPath(SeqeraFileSystem fs, String uriString) {
        this.fs = fs
        this.relPath = null
        if (!uriString.startsWith("${SCHEME}://"))
            throw new InvalidPathException(uriString, "Not a seqera:// URI")
        // strip scheme: seqera://rest
        final withoutScheme = uriString.substring("${SCHEME}://".length())
        // split on '/'
        final parts = withoutScheme.split('/', -1).toList().findAll { it != null } as List<String>
        // parts[0]=org, parts[1]=workspace, parts[2]=resourceType, parts[3]=datasetName[@version]
        this.org = parts.size() > 0 && parts[0] ? parts[0] : null
        this.workspace = parts.size() > 1 && parts[1] ? parts[1] : null
        this.resourceType = parts.size() > 2 && parts[2] ? parts[2] : null
        if (parts.size() > 3 && parts[3]) {
            final last = parts[3]
            final atIdx = last.lastIndexOf('@')
            if (atIdx > 0) {
                this.datasetName = last.substring(0, atIdx)
                this.version = last.substring(atIdx + 1)
            } else {
                this.datasetName = last
                this.version = null
            }
        } else {
            this.datasetName = null
            this.version = null
        }
        validatePath(uriString)
    }

    /** Internal constructor for programmatic absolute path creation */
    SeqeraPath(SeqeraFileSystem fs, String org, String workspace, String resourceType, String datasetName, String version) {
        this.fs = fs
        this.relPath = null
        this.org = org
        this.workspace = workspace
        this.resourceType = resourceType
        this.datasetName = datasetName
        this.version = version
        validatePath(null)
    }

    /**
     * Constructor for relative paths produced by {@link #relativize(Path)}.
     * The {@code relPath} is a slash-separated string of the differing path segments.
     * All segment fields are {@code null}; {@link #isAbsolute()} returns {@code false}.
     */
    SeqeraPath(String relPath) {
        this.fs = null
        this.relPath = relPath ?: ''
        this.org = null
        this.workspace = null
        this.resourceType = null
        this.datasetName = null
        this.version = null
    }

    /**
     * Validate structural integrity: deeper segments require all shallower ones,
     * and no segment may contain {@code /}.
     *
     * @param original original URI string used in error messages (null → derive from fields)
     * @throws InvalidPathException if the path is malformed
     */
    private void validatePath(String original) {
        final label = original ?: rawPath()
        if (datasetName && !workspace)
            throw new InvalidPathException(label, "Dataset path requires a workspace segment")
        if (resourceType && !workspace)
            throw new InvalidPathException(label, "Resource type requires a workspace segment")
        if (workspace && !org)
            throw new InvalidPathException(label, "Workspace requires an org segment")
        // Segments from URI parsing never contain '/', but guard the internal constructor too
        if (org?.contains('/'))
            throw new InvalidPathException(label, "Org name cannot contain '/'")
        if (workspace?.contains('/'))
            throw new InvalidPathException(label, "Workspace name cannot contain '/'")
        if (resourceType?.contains('/'))
            throw new InvalidPathException(label, "Resource type cannot contain '/'")
        if (datasetName?.contains('/'))
            throw new InvalidPathException(label, "Dataset name cannot contain '/'")
    }

    /** Return a list of name component strings (works for both absolute and relative paths). */
    private List<String> nameComponents() {
        if (isAbsolute()) {
            final d = depth()
            final result = new ArrayList<String>(d)
            for (int i = 0; i < d; i++)
                result.add(getName(i).toString())
            return result
        }
        if (!relPath) return Collections.<String>emptyList()
        return relPath.split('/').toList().findAll { String s -> s } as List<String>
    }

    /** Build a raw path string from the current fields, for use in exception messages. */
    private String rawPath() {
        final sb = new StringBuilder("${SCHEME}://")
        if (org) sb.append(org)
        if (workspace) sb.append('/').append(workspace)
        if (resourceType) sb.append('/').append(resourceType)
        if (datasetName) {
            sb.append('/').append(datasetName)
            if (version) sb.append('@').append(version)
        }
        return sb.toString()
    }

    // ---- path component accessors ----

    String getOrg() { org }
    String getWorkspace() { workspace }
    String getResourceType() { resourceType }
    String getDatasetName() { datasetName }
    String getVersion() { version }

    /**
     * Path depth: 0=root, 1=org, 2=workspace, 3=resourceType, 4=dataset file.
     */
    int depth() {
        if (datasetName) return 4
        if (resourceType) return 3
        if (workspace) return 2
        if (org) return 1
        return 0
    }

    boolean isDirectory() { depth() < 4 }
    boolean isRegularFile() { depth() == 4 }

    // ---- Path API ----

    @Override
    FileSystem getFileSystem() { fs }

    @Override
    boolean isAbsolute() { fs != null }

    @Override
    Path getRoot() { new SeqeraPath(fs, null, null, null, null, null) }

    @Override
    Path getFileName() {
        final d = depth()
        if (d == 0) return null
        final name = d == 4 ? (version ? "${datasetName}@${version}" : datasetName)
                   : d == 3 ? resourceType
                   : d == 2 ? workspace
                   : org
        return new SeqeraPath( name as String)
    }

    @Override
    Path getParent() {
        final d = depth()
        if (d == 0) return null
        if (d == 1) return new SeqeraPath(fs, null, null, null, null, null)
        if (d == 2) return new SeqeraPath(fs, org, null, null, null, null)
        if (d == 3) return new SeqeraPath(fs, org, workspace, null, null, null)
        return new SeqeraPath(fs, org, workspace, resourceType, null, null)
    }

    @Override
    int getNameCount() { depth() }

    @Override
    Path getName(int index) {
        final d = depth()
        if (index < 0 || index >= d)
            throw new IllegalArgumentException("Index out of range: $index")
        if (index == 0) return new SeqeraPath(org)
        if (index == 1) return new SeqeraPath(workspace)
        if (index == 2) return new SeqeraPath(resourceType)
        return new SeqeraPath((version ? "${datasetName}@${version}" : datasetName) as String)
    }

    @Override
    Path subpath(int beginIndex, int endIndex) {
        throw new UnsupportedOperationException("subpath not supported by seqera:// paths")
    }

    @Override
    boolean startsWith(Path other) {
        if (other !instanceof SeqeraPath)
            return false
        final that = (SeqeraPath) other
        if (this.isAbsolute() != that.isAbsolute())
            return false
        final thisNames = this.nameComponents()
        final thatNames = that.nameComponents()
        if (thatNames.size() > thisNames.size())
            return false
        for (int i = 0; i < thatNames.size(); i++) {
            if (thisNames[i] != thatNames[i])
                return false
        }
        return true
    }

    @Override
    boolean startsWith(String other) {
        if (!other) return false
        try {
            final Path p = SeqeraPath.isSeqeraUri(other) ? new SeqeraPath(fs, other) : new SeqeraPath(other)
            return startsWith(p)
        } catch (Exception ignored) {
            return false
        }
    }

    @Override
    boolean endsWith(Path other) {
        if (other !instanceof SeqeraPath)
            return false
        final that = (SeqeraPath) other
        if (that.isAbsolute())
            return this.equals(that)
        final thisNames = this.nameComponents()
        final thatNames = that.nameComponents()
        if (thatNames.isEmpty() || thatNames.size() > thisNames.size())
            return false
        final offset = thisNames.size() - thatNames.size()
        for (int i = 0; i < thatNames.size(); i++) {
            if (thisNames[offset + i] != thatNames[i])
                return false
        }
        return true
    }

    @Override
    boolean endsWith(String other) {
        if (!other) return false
        try {
            final Path p = SeqeraPath.isSeqeraUri(other) ? new SeqeraPath(fs, other) : new SeqeraPath(other)
            return endsWith(p)
        } catch (Exception ignored) {
            return false
        }
    }

    @Override
    Path normalize() { this }

    @Override
    Path resolve(Path other) {
        if (other instanceof SeqeraPath) {
            final that = (SeqeraPath) other
            if (that.isAbsolute()) return that
            // Relative SeqeraPath: resolve each segment of relPath against this
            return resolve(that.relPath)
        }
        return resolve(other.toString())
    }

    @Override
    Path resolve(String segment) {
        if (!segment) return this
        // Absolute seqera:// URI — parse and return directly
        if (segment.startsWith(PROTOCOL))
            return new SeqeraPath(fs, segment)
        // Strip a single leading slash
        final stripped = segment.startsWith(SEPARATOR) ? segment.substring(1) : segment
        if (!stripped) return this
        // Multi-segment: split and resolve one segment at a time
        final segs = stripped.split(SEPARATOR, -1).findAll { String s -> s } as List<String>
        SeqeraPath result = this
        for (String seg : segs) {
            result = result.resolveOne(seg)
        }
        return result
    }

    /** Resolve a single (non-empty, slash-free) segment against this path. */
    private SeqeraPath resolveOne(String seg) {
        final d = depth()
        if (d == 0) return new SeqeraPath(fs, seg, null, null, null, null)
        if (d == 1) return new SeqeraPath(fs, org, seg, null, null, null)
        if (d == 2) return new SeqeraPath(fs, org, workspace, seg, null, null)
        if (d == 3) {
            final atIdx = seg.lastIndexOf('@')
            if (atIdx > 0)
                return new SeqeraPath(fs, org, workspace, resourceType, seg.substring(0, atIdx), seg.substring(atIdx + 1))
            return new SeqeraPath(fs, org, workspace, resourceType, seg, null)
        }
        throw new IllegalStateException("Cannot resolve a path segment on a depth-4 path: $this")
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
        if (other !instanceof SeqeraPath)
            throw new ProviderMismatchException()
        final that = (SeqeraPath) other
        if (!this.isAbsolute() || !that.isAbsolute())
            throw new IllegalArgumentException("Both paths must be absolute to relativize: ${this} vs ${other}")
        final thisNames = this.nameComponents()
        final thatNames = that.nameComponents()
        // Find common prefix length
        int common = 0
        while (common < thisNames.size() && common < thatNames.size()
               && thisNames[common] == thatNames[common])
            common++
        // Build ".." for each remaining segment in this, then append remaining segments of other
        final parts = new ArrayList<String>()
        for (int i = common; i < thisNames.size(); i++)
            parts.add('..')
        for (int i = common; i < thatNames.size(); i++)
            parts.add(thatNames[i])
        return new SeqeraPath(parts.join(SEPARATOR))
    }

    @Override
    URI toUri() {
        // Build path component for depth >= 2
        String uriPath = null
        if (workspace) {
            final segments = [workspace]
            if (resourceType) segments.add(resourceType)
            if (datasetName) segments.add(version ? "${datasetName}@${version}" as String : datasetName)
            uriPath = '/' + segments.join('/')
        }
        // new URI(scheme, authority, path, query, fragment) avoids URI.create() pitfalls for edge cases
        return new URI(SCHEME, org ?: '', uriPath, null, null)
    }

    @Override
    String toString() {
        if (!isAbsolute()) return relPath
        // Return the canonical human-readable representation
        final d = depth()
        if (d == 0) return "${SCHEME}://"
        return toUri().toString()
    }

    @Override
    Path toAbsolutePath() {
        if (!isAbsolute())
            throw new IllegalStateException("Cannot convert relative SeqeraPath to absolute — no default directory context")
        return this
    }

    @Override
    Path toRealPath(LinkOption... options) { this }

    @Override
    File toFile() {
        throw new UnsupportedOperationException("toFile() not supported for seqera:// paths")
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>[] events, WatchEvent.Modifier... modifiers) {
        throw new UnsupportedOperationException("WatchService not supported by seqera:// paths")
    }

    @Override
    WatchKey register(WatchService watcher, WatchEvent.Kind<?>... events) {
        throw new UnsupportedOperationException("WatchService not supported by seqera:// paths")
    }

    @Override
    Iterator<Path> iterator() {
        final d = depth()
        final List<Path> parts = new ArrayList<>(d)
        for (int i = 0; i < d; i++) {
            parts.add(getName(i))
        }
        return parts.iterator()
    }

    @Override
    int compareTo(Path other) {
        return toString().compareTo(other.toString())
    }

    @Override
    boolean equals(Object obj) {
        if (obj == this) return true
        if (obj !instanceof SeqeraPath) return false
        return toString() == obj.toString()
    }

    @Override
    int hashCode() { toString().hashCode() }

    static URI asUri(String path) {
        if( !path )
            throw new IllegalArgumentException("Missing 'path' argument")
        if( !path.startsWith(PROTOCOL) )
            throw new IllegalArgumentException("Invalid Seqera file system path URI - it must start with '${PROTOCOL}' prefix - offending value: $path")
        if( path.startsWith(PROTOCOL + SEPARATOR) && path.length() > PROTOCOL.length() + 1 )
            throw new IllegalArgumentException("Invalid Seqera file system path URI - make sure the scheme prefix does not contain more than two slash characters or a query in the root '/' - offending value: $path")

        //URI strings like seqera://./something are converted to seqera://something
        if( path.startsWith(PROTOCOL + './') ) {
            path = PROTOCOL + path.substring(PROTOCOL.length() + 2)
        }

        if( path == PROTOCOL || path == PROTOCOL + '.') //Empty path case
            return new URI(PROTOCOL + '/')
        return new URI(path)
    }

    static boolean isSeqeraUri(String path) {
        return path && path.startsWith(PROTOCOL)
    }
}
