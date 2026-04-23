# Tasks: Seqera NIO Filesystem Support for Platform Data-Links

**Branch**: `260422-seqera-datalinks-fs` | **Spec**: [spec.md](spec.md) | **Plan**: [plan.md](plan.md)

> **For agentic workers**: execute tasks in order. Each task is self-contained and ends with a commit step. Do not skip TDD steps — write the test first, watch it fail, then make it pass. All commits use `git commit -s`.

Tests use Spock with `Mock(TowerClient)` + `groovy.json.JsonOutput` fixtures — matching the style of `SeqeraDatasetClientTest` and `SeqeraFileSystemProviderTest`. No WireMock, no real HTTP.

Legend:
- **[P]**: can be done in parallel with the previous task (different files, no dependency)
- **[Story]**: which user story from the spec
- Exact file paths are relative to repo root

---

## Phase 1: Foundational Refactor (blocks all US)

**Purpose**: Extract `ResourceTypeHandler`, move dataset-specific logic out of the generic classes. Existing dataset tests must pass unchanged at the end of this phase.

### T001 — Generalize `SeqeraFileAttributes`

**Files:**
- Modify: `plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/SeqeraFileAttributes.groovy`

- [ ] **Step 1: Replace the DatasetDto-coupled class with a generic one**

Replace the entire class body (everything between `class SeqeraFileAttributes implements BasicFileAttributes {` and the final `}`) with:

```groovy
@CompileStatic
class SeqeraFileAttributes implements BasicFileAttributes {

    private final boolean directory
    private final long size
    private final Instant lastModified
    private final Instant created
    private final Object fileKey

    /** Construct attributes for a virtual directory (any depth). */
    SeqeraFileAttributes(boolean isDir) {
        this.directory = isDir
        this.size = 0L
        this.lastModified = Instant.EPOCH
        this.created = Instant.EPOCH
        this.fileKey = null
    }

    /** Construct attributes for a regular file. */
    SeqeraFileAttributes(long size, Instant lastModified, Instant created, Object fileKey) {
        this.directory = false
        this.size = size ?: 0L
        this.lastModified = lastModified ?: Instant.EPOCH
        this.created = created ?: Instant.EPOCH
        this.fileKey = fileKey
    }

    @Override FileTime lastModifiedTime() { FileTime.from(lastModified) }
    @Override FileTime lastAccessTime() { FileTime.from(lastModified) }
    @Override FileTime creationTime() { FileTime.from(created) }
    @Override boolean isRegularFile() { !directory }
    @Override boolean isDirectory() { directory }
    @Override boolean isSymbolicLink() { false }
    @Override boolean isOther() { false }
    @Override long size() { size }
    @Override Object fileKey() { fileKey }
}
```

- [ ] **Step 2: Drop the now-unused `DatasetDto` import**

In the same file, remove `import io.seqera.tower.model.DatasetDto`. Keep `import java.time.Instant`.

- [ ] **Step 3: Compile the plugin**

Run: `./gradlew :plugins:nf-tower:compileGroovy`
Expected: `BUILD FAILED` with errors in `SeqeraFileSystemProvider.groovy` (calls `new SeqeraFileAttributes(dataset)`) — we will fix that when we extract the dataset handler. Compile errors here are expected at this step only; leave them for T004.

- [ ] **Step 4: Commit**

```bash
git add plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/SeqeraFileAttributes.groovy
git commit -s -m "refactor(nf-tower): generalize SeqeraFileAttributes (no DatasetDto coupling)"
```

### T002 — Add `ResourceTypeHandler` interface

**Files:**
- Create: `plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/ResourceTypeHandler.groovy`

- [ ] **Step 1: Write the interface**

```groovy
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

import java.nio.file.AccessMode
import java.nio.file.Path

/**
 * Strategy owning the semantics of one depth-3 path segment under {@code seqera://}.
 * Registered in {@link SeqeraFileSystem} at filesystem construction.
 */
interface ResourceTypeHandler {

    /** e.g. {@code "datasets"} or {@code "data-links"}. Must match the depth-3 path segment. */
    String getResourceType()

    /** List entries at the given directory path. Caller has verified depth &ge; 3 and {@code sp.isDirectory()}. */
    List<Path> list(SeqeraPath dir) throws IOException

    /** Return attributes for any path at depth &ge; 3 owned by this handler. */
    SeqeraFileAttributes readAttributes(SeqeraPath path) throws IOException

    /** Open a read stream for a leaf path. Throw {@link IllegalArgumentException} if the path is a directory. */
    InputStream newInputStream(SeqeraPath path) throws IOException

    /** Verify the path exists and requested modes are satisfiable. READ is allowed; WRITE/EXECUTE throw {@link java.nio.file.AccessDeniedException}. */
    void checkAccess(SeqeraPath path, AccessMode... modes) throws IOException
}
```

- [ ] **Step 2: Compile (the other existing errors from T001 still stand — expected)**

Run: `./gradlew :plugins:nf-tower:compileGroovy`
Expected: same errors as T001 (in `SeqeraFileSystemProvider.groovy`), no new errors from the new file.

- [ ] **Step 3: Commit**

```bash
git add plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/ResourceTypeHandler.groovy
git commit -s -m "refactor(nf-tower): add ResourceTypeHandler interface"
```

### T003 — Refactor `SeqeraPath` to use generic `trail` segments

**Files:**
- Modify: `plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/SeqeraPath.groovy`
- Modify: `plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/SeqeraPathTest.groovy`

- [ ] **Step 1: Read the current `SeqeraPath` tests to understand coverage**

Open `plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/SeqeraPathTest.groovy`. Note the existing cases you must continue to pass (URI parsing, `toUri()` round-trip, `getParent`, `resolve`, `@version`, `equals`/`hashCode`).

- [ ] **Step 2: Write failing tests for the new depth-≥5 cases**

Append to the end of the class body (before the final `}`):

```groovy
    def "parse data-link path with provider and name"() {
        when:
        def p = new SeqeraPath(Mock(SeqeraFileSystem), 'seqera://acme/research/data-links/aws/inputs')

        then:
        p.org == 'acme'
        p.workspace == 'research'
        p.resourceType == 'data-links'
        p.trail == ['aws', 'inputs']
        p.depth() == 5
    }

    def "parse data-link path with nested sub-path"() {
        when:
        def p = new SeqeraPath(Mock(SeqeraFileSystem), 'seqera://acme/research/data-links/aws/inputs/reads/sample1.fq.gz')

        then:
        p.trail == ['aws', 'inputs', 'reads', 'sample1.fq.gz']
        p.depth() == 7
        p.isRegularFile() == false  // handler decides; generic class reports directory-by-default for trail != 1 dataset
    }

    def "getParent walks up one trail segment for deep data-link paths"() {
        given:
        def fs = Mock(SeqeraFileSystem)
        def p = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs/reads/s.fq.gz')

        when:
        def parent = p.getParent() as SeqeraPath

        then:
        parent.trail == ['aws', 'inputs', 'reads']
        parent.depth() == 6
    }

    def "resolve appends one segment to trail"() {
        given:
        def fs = Mock(SeqeraFileSystem)
        def p = new SeqeraPath(fs, 'seqera://acme/research/data-links/aws/inputs')

        when:
        def child = p.resolve('reads') as SeqeraPath

        then:
        child.trail == ['aws', 'inputs', 'reads']
    }

    def "toUri round-trip for deep data-link path"() {
        given:
        def fs = Mock(SeqeraFileSystem)
        def original = 'seqera://acme/research/data-links/aws/inputs/reads/sample.fq.gz'

        when:
        def p = new SeqeraPath(fs, original)

        then:
        p.toUri().toString() == original
    }

    def "dataset version pinning preserved after refactor"() {
        given:
        def fs = Mock(SeqeraFileSystem)

        when:
        def p = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples@3')

        then:
        p.resourceType == 'datasets'
        p.trail == ['samples']
        p.version == '3'
        p.datasetName == 'samples'
        p.depth() == 4
    }
```

Run: `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.SeqeraPathTest' -i`
Expected: these new tests fail (methods `trail`, `version`, `datasetName` don't behave right yet, or `depth()` returns 4 for data-link paths because the current parse truncates).

- [ ] **Step 3: Rewrite `SeqeraPath` to use generic trail segments**

Replace `SeqeraPath.groovy` entirely with:

```groovy
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
 * The generic class is resource-type-agnostic for depth &ge; 3: segments after
 * {@code resourceType} are exposed as {@link #getTrail()} for the matching
 * {@link ResourceTypeHandler} to interpret.
 *
 * The dataset convention (single trail segment, optional {@code @version} suffix)
 * is preserved via {@link #getDatasetName()} and {@link #getVersion()} accessors.
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
    private final String version
    private final String relPath

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
        final List<String> tail = parts.size() > 3 ? new ArrayList<String>(parts.subList(3, parts.size())) : new ArrayList<String>()
        // For datasets: strip "@version" from the last trail segment if present.
        if (this.resourceType == 'datasets' && tail.size() == 1) {
            final last = tail[0]
            final atIdx = last.lastIndexOf('@')
            if (atIdx > 0) {
                tail[0] = last.substring(0, atIdx)
                this.version = last.substring(atIdx + 1)
            } else {
                this.version = null
            }
        } else {
            this.version = null
        }
        this.trail = Collections.unmodifiableList(tail)
        validatePath(uriString)
    }

    /** Programmatic absolute-path constructor. */
    SeqeraPath(SeqeraFileSystem fs, String org, String workspace, String resourceType, List<String> trail, String version) {
        this.fs = fs
        this.relPath = null
        this.org = org
        this.workspace = workspace
        this.resourceType = resourceType
        this.trail = trail != null ? Collections.unmodifiableList(new ArrayList<String>(trail)) : Collections.<String>emptyList()
        this.version = version
        validatePath(null)
    }

    /** Relative path, produced only by {@link #relativize(Path)}. */
    SeqeraPath(String relPath) {
        this.fs = null
        this.relPath = relPath ?: ''
        this.org = null
        this.workspace = null
        this.resourceType = null
        this.trail = Collections.<String>emptyList()
        this.version = null
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
        // Datasets accept at most one trail segment
        if (resourceType == 'datasets' && trail.size() > 1)
            throw new InvalidPathException(label, "Dataset paths cannot have sub-paths beyond the dataset name")
    }

    private String rawPath() {
        final sb = new StringBuilder(PROTOCOL)
        if (org) sb.append(org)
        if (workspace) sb.append('/').append(workspace)
        if (resourceType) sb.append('/').append(resourceType)
        for (int i = 0; i < trail.size(); i++) {
            sb.append('/')
            if (i == trail.size() - 1 && version)
                sb.append(trail[i]).append('@').append(version)
            else
                sb.append(trail[i])
        }
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
    String getVersion() { version }

    /** Backwards-compat: dataset name is the single trail segment when resourceType=='datasets'. */
    String getDatasetName() {
        (resourceType == 'datasets' && trail.size() == 1) ? trail[0] : null
    }

    int depth() {
        if (resourceType) return 3 + trail.size()
        if (workspace) return 2
        if (org) return 1
        return 0
    }

    boolean isDirectory() {
        // Dataset leaf at depth 4 is a file; all other shapes are directory-by-default.
        // Handlers override this interpretation for data-link sub-paths via readAttributes.
        !(resourceType == 'datasets' && trail.size() == 1)
    }

    boolean isRegularFile() { !isDirectory() }

    // ---- Path API ----

    @Override FileSystem getFileSystem() { fs }
    @Override boolean isAbsolute() { fs != null }

    @Override
    Path getRoot() { new SeqeraPath(fs, null, null, null, null, null) }

    @Override
    Path getFileName() {
        final d = depth()
        if (d == 0) return null
        if (d >= 4) {
            final last = trail[trail.size() - 1]
            return new SeqeraPath((d == 4 && version) ? "${last}@${version}" as String : last)
        }
        if (d == 3) return new SeqeraPath(resourceType)
        if (d == 2) return new SeqeraPath(workspace)
        return new SeqeraPath(org)
    }

    @Override
    Path getParent() {
        final d = depth()
        if (d == 0) return null
        if (d == 1) return new SeqeraPath(fs, null, null, null, null, null)
        if (d == 2) return new SeqeraPath(fs, org, null, null, null, null)
        if (d == 3) return new SeqeraPath(fs, org, workspace, null, null, null)
        // d >= 4: drop last trail segment
        final newTrail = trail.subList(0, trail.size() - 1)
        return new SeqeraPath(fs, org, workspace, resourceType, newTrail, null)
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
        final trailIdx = index - 3
        final seg = trail[trailIdx]
        // Only the last segment of a depth-4 dataset path carries the version suffix
        if (trailIdx == trail.size() - 1 && version && resourceType == 'datasets')
            return new SeqeraPath("${seg}@${version}" as String)
        return new SeqeraPath(seg)
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
            final p = SeqeraPath.isSeqeraUri(other) ? new SeqeraPath(fs, other) : new SeqeraPath(other)
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
            final p = SeqeraPath.isSeqeraUri(other) ? new SeqeraPath(fs, other) : new SeqeraPath(other)
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
        if (d == 0) return new SeqeraPath(fs, seg, null, null, null, null)
        if (d == 1) return new SeqeraPath(fs, org, seg, null, null, null)
        if (d == 2) return new SeqeraPath(fs, org, workspace, seg, null, null)
        // d >= 3: append to trail (with @version parsing only for dataset-shaped paths at depth 3→4)
        if (d == 3 && resourceType == 'datasets') {
            final atIdx = seg.lastIndexOf('@')
            if (atIdx > 0)
                return new SeqeraPath(fs, org, workspace, resourceType, [seg.substring(0, atIdx)], seg.substring(atIdx + 1))
        }
        final newTrail = new ArrayList<String>(trail)
        newTrail.add(seg)
        return new SeqeraPath(fs, org, workspace, resourceType, newTrail, null)
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
            for (int i = 0; i < trail.size(); i++) {
                final t = trail[i]
                if (i == trail.size() - 1 && version && resourceType == 'datasets')
                    segments.add("${t}@${version}" as String)
                else
                    segments.add(t)
            }
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
```

- [ ] **Step 4: Run the `SeqeraPathTest` — new and existing cases must pass**

Run: `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.SeqeraPathTest' -i`
Expected: all tests pass. If any existing test fails, it signals a refactor regression — fix before continuing.

- [ ] **Step 5: Commit**

```bash
git add plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/SeqeraPath.groovy \
        plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/SeqeraPathTest.groovy
git commit -s -m "refactor(nf-tower): generalize SeqeraPath with trail segments for multi-resource support"
```

### T004 — Extract `DatasetsResourceHandler`

**Files:**
- Create: `plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/handler/DatasetsResourceHandler.groovy`
- Create: `plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/handler/DatasetsResourceHandlerTest.groovy`

- [ ] **Step 1: Write the handler — move logic from current `SeqeraFileSystemProvider`**

The source of the logic is the current `newInputStream`, `readAttributes`, `newDirectoryStream` (depth 3 only), and `checkAccess` branches in `SeqeraFileSystemProvider.groovy` that handle `datasets`. Preserve behavior byte-for-byte.

```groovy
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

package io.seqera.tower.plugin.fs.handler

import java.nio.file.AccessDeniedException
import java.nio.file.AccessMode
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.time.Instant

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.model.DatasetDto
import io.seqera.tower.model.DatasetVersionDto
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient
import io.seqera.tower.plugin.fs.ResourceTypeHandler
import io.seqera.tower.plugin.fs.SeqeraFileAttributes
import io.seqera.tower.plugin.fs.SeqeraFileSystem
import io.seqera.tower.plugin.fs.SeqeraPath

/**
 * {@link ResourceTypeHandler} for {@code datasets} resource type.
 * All logic previously inlined in {@link SeqeraFileSystemProvider} for dataset paths lives here.
 */
@Slf4j
@CompileStatic
class DatasetsResourceHandler implements ResourceTypeHandler {

    public static final String TYPE = 'datasets'

    private final SeqeraFileSystem fs
    private final SeqeraDatasetClient client

    DatasetsResourceHandler(SeqeraFileSystem fs, SeqeraDatasetClient client) {
        this.fs = fs
        this.client = client
    }

    @Override
    String getResourceType() { TYPE }

    @Override
    List<Path> list(SeqeraPath dir) throws IOException {
        final d = dir.depth()
        if (d == 3) {
            final workspaceId = fs.resolveWorkspaceId(dir.org, dir.workspace)
            final datasets = fs.resolveDatasets(workspaceId)
            return datasets.collect { DatasetDto ds -> dir.resolve(ds.name) as Path }
        }
        throw new IllegalArgumentException("datasets handler cannot list depth $d paths: $dir")
    }

    @Override
    SeqeraFileAttributes readAttributes(SeqeraPath p) throws IOException {
        final d = p.depth()
        if (d == 3) {
            // resource-type dir — validate workspace
            fs.resolveWorkspaceId(p.org, p.workspace)
            return new SeqeraFileAttributes(true)
        }
        if (d != 4)
            throw new NoSuchFileException(p.toString(), null, "Invalid dataset path depth: $d")
        final workspaceId = fs.resolveWorkspaceId(p.org, p.workspace)
        final dataset = fs.resolveDataset(workspaceId, p.datasetName)
        if (!dataset)
            throw new NoSuchFileException(p.toString(), null, "Dataset '${p.datasetName}' not found in workspace ${p.workspace}")
        return new SeqeraFileAttributes(
            0L,
            dataset.lastUpdated?.toInstant() ?: Instant.EPOCH,
            dataset.dateCreated?.toInstant() ?: Instant.EPOCH,
            dataset.id )
    }

    @Override
    InputStream newInputStream(SeqeraPath p) throws IOException {
        if (p.depth() != 4)
            throw new IllegalArgumentException("Operation `newInputStream` requires a dataset path (depth 4): $p")
        final workspaceId = fs.resolveWorkspaceId(p.org, p.workspace)
        final dataset = fs.resolveDataset(workspaceId, p.datasetName)
        if (!dataset)
            throw new NoSuchFileException(p.toString(), null, "Dataset '${p.datasetName}' not found in workspace ${p.workspace}")
        final version = resolveVersion(dataset, p)
        log.debug "Downloading dataset '${p.datasetName}' version ${version.version} (${version.fileName}) from workspace $workspaceId"
        return client.downloadDataset(dataset.id, String.valueOf(version.version), version.fileName, dataset.workspaceId)
    }

    @Override
    void checkAccess(SeqeraPath p, AccessMode... modes) throws IOException {
        for (AccessMode m : modes) {
            if (m == AccessMode.WRITE || m == AccessMode.EXECUTE)
                throw new AccessDeniedException(p.toString(), null, "seqera:// datasets are read-only")
        }
        // READ: make sure the dataset resolves
        readAttributes(p)
    }

    private DatasetVersionDto resolveVersion(DatasetDto dataset, SeqeraPath p) throws IOException {
        final pinned = p.version
        final versions = fs.resolveVersions(dataset.id, dataset.workspaceId)
        if (versions.isEmpty())
            throw new NoSuchFileException(p.toString(), null, "No versions available for dataset '${dataset.name}'")
        if (pinned) {
            final found = versions.find { DatasetVersionDto v -> String.valueOf(v.version) == pinned }
            if (!found)
                throw new NoSuchFileException(p.toString(), null, "Version '$pinned' not found for dataset '${dataset.name}'")
            return found
        }
        final latest = versions.findAll { DatasetVersionDto v -> !v.disabled }.max { DatasetVersionDto v -> v.version }
        if (!latest)
            throw new NoSuchFileException(p.toString(), null, "No enabled versions for dataset '${dataset.name}'")
        return latest
    }
}
```

- [ ] **Step 2: Write Spock tests for the handler**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * (Apache 2.0)
 */
package io.seqera.tower.plugin.fs.handler

import java.nio.file.AccessMode
import java.nio.file.NoSuchFileException

import io.seqera.tower.model.DatasetDto
import io.seqera.tower.model.DatasetVersionDto
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient
import io.seqera.tower.plugin.fs.SeqeraFileSystem
import io.seqera.tower.plugin.fs.SeqeraPath
import spock.lang.Specification

class DatasetsResourceHandlerTest extends Specification {

    def fs = Mock(SeqeraFileSystem)
    def client = Mock(SeqeraDatasetClient)
    def handler = new DatasetsResourceHandler(fs, client)

    def "getResourceType returns 'datasets'"() {
        expect:
        handler.resourceType == 'datasets'
    }

    def "list at depth 3 returns one path per dataset"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets')
        def ds1 = new DatasetDto(id: 'd1', name: 'one')
        def ds2 = new DatasetDto(id: 'd2', name: 'two')

        when:
        def paths = handler.list(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * fs.resolveDatasets(10L) >> [ds1, ds2]
        paths*.toString() == ['seqera://acme/research/datasets/one', 'seqera://acme/research/datasets/two']
    }

    def "newInputStream resolves latest non-disabled version when no pin"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')
        def ds = new DatasetDto(id: 'd1', name: 'samples', workspaceId: 10L)
        def v1 = new DatasetVersionDto(datasetId: 'd1', version: 1L, fileName: 'a.csv', disabled: false)
        def v2 = new DatasetVersionDto(datasetId: 'd1', version: 2L, fileName: 'b.csv', disabled: false)
        def content = new ByteArrayInputStream('x'.bytes)

        when:
        def stream = handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * fs.resolveDataset(10L, 'samples') >> ds
        1 * fs.resolveVersions('d1', 10L) >> [v1, v2]
        1 * client.downloadDataset('d1', '2', 'b.csv', 10L) >> content
        stream === content
    }

    def "newInputStream honors @version pin"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples@1')
        def ds = new DatasetDto(id: 'd1', name: 'samples', workspaceId: 10L)
        def v1 = new DatasetVersionDto(datasetId: 'd1', version: 1L, fileName: 'a.csv', disabled: false)
        def v2 = new DatasetVersionDto(datasetId: 'd1', version: 2L, fileName: 'b.csv', disabled: false)

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * fs.resolveDataset(10L, 'samples') >> ds
        1 * fs.resolveVersions('d1', 10L) >> [v1, v2]
        1 * client.downloadDataset('d1', '1', 'a.csv', 10L) >> new ByteArrayInputStream('x'.bytes)
    }

    def "newInputStream throws NoSuchFileException when dataset is missing"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/ghost')

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * fs.resolveDataset(10L, 'ghost') >> null
        thrown(NoSuchFileException)
    }

    def "checkAccess rejects WRITE"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/datasets/samples')

        when:
        handler.checkAccess(path, AccessMode.WRITE)

        then:
        thrown(java.nio.file.AccessDeniedException)
    }
}
```

- [ ] **Step 3: Run the handler tests**

Run: `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.handler.DatasetsResourceHandlerTest' -i`
Expected: all pass (after T005–T006 make the provider compile; until then the top-level compile failure blocks this — continue to T005).

- [ ] **Step 4: Commit**

```bash
git add plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/handler/DatasetsResourceHandler.groovy \
        plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/handler/DatasetsResourceHandlerTest.groovy
git commit -s -m "refactor(nf-tower): extract DatasetsResourceHandler from SeqeraFileSystemProvider"
```

### T005 — Add handler registry to `SeqeraFileSystem`

**Files:**
- Modify: `plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/SeqeraFileSystem.groovy`

- [ ] **Step 1: Add the registry field and accessors**

Insert after the existing cache field declarations (e.g. after the `versionCache` line):

```groovy
    private final Map<String, ResourceTypeHandler> handlers = new LinkedHashMap<>()
```

Then insert immediately before the closing `}` of the class:

```groovy
    // ---- handler registry ----

    synchronized void registerHandler(ResourceTypeHandler handler) {
        handlers.put(handler.resourceType, handler)
    }

    synchronized ResourceTypeHandler getHandler(String resourceType) {
        handlers.get(resourceType)
    }

    synchronized Set<String> getResourceTypes() {
        Collections.unmodifiableSet(new LinkedHashSet<String>(handlers.keySet()))
    }
```

- [ ] **Step 2: Compile**

Run: `./gradlew :plugins:nf-tower:compileGroovy`
Expected: compile errors in `SeqeraFileSystemProvider.groovy` remain; no new errors from `SeqeraFileSystem`.

- [ ] **Step 3: Commit**

```bash
git add plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/SeqeraFileSystem.groovy
git commit -s -m "refactor(nf-tower): add ResourceTypeHandler registry to SeqeraFileSystem"
```

### T006 — Refactor `SeqeraFileSystemProvider` to dispatch via handlers

**Files:**
- Modify: `plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/SeqeraFileSystemProvider.groovy`

- [ ] **Step 1: Replace the provider body**

Replace the class body with the following (the class shell, package, imports are adjusted):

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * (Apache 2.0 header unchanged)
 */

package io.seqera.tower.plugin.fs

import java.nio.channels.SeekableByteChannel
import java.nio.file.AccessDeniedException
import java.nio.file.AccessMode
import java.nio.file.CopyOption
import java.nio.file.DirectoryIteratorException
import java.nio.file.DirectoryStream
import java.nio.file.FileStore
import java.nio.file.FileSystem
import java.nio.file.FileSystemAlreadyExistsException
import java.nio.file.FileSystemNotFoundException
import java.nio.file.Files
import java.nio.file.LinkOption
import java.nio.file.NoSuchFileException
import java.nio.file.NotDirectoryException
import java.nio.file.OpenOption
import java.nio.file.Path
import java.nio.file.ProviderMismatchException
import java.nio.file.StandardOpenOption
import java.nio.file.attribute.BasicFileAttributes
import java.nio.file.attribute.FileAttribute
import java.nio.file.attribute.FileAttributeView
import java.nio.file.spi.FileSystemProvider

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.plugin.TowerClient
import io.seqera.tower.plugin.TowerFactory
import io.seqera.tower.plugin.dataset.SeqeraDatasetClient
import io.seqera.tower.plugin.datalink.SeqeraDataLinkClient
import io.seqera.tower.plugin.fs.handler.DatasetsResourceHandler
import io.seqera.tower.plugin.fs.handler.DataLinksResourceHandler

@Slf4j
@CompileStatic
class SeqeraFileSystemProvider extends FileSystemProvider {

    public static final String SCHEME = 'seqera'

    private volatile SeqeraFileSystem fileSystem

    @Override String getScheme() { SCHEME }

    @Override
    synchronized FileSystem newFileSystem(URI uri, Map<String, ?> env) throws IOException {
        checkScheme(uri)
        if (fileSystem)
            throw new FileSystemAlreadyExistsException("File system `seqera://` already exists")
        final TowerClient tc = TowerFactory.client()
        if (!tc)
            throw new IllegalStateException("File system `seqera://` requires the Seqera Platform access token — use `tower.accessToken` or TOWER_ACCESS_TOKEN")
        final datasetClient = new SeqeraDatasetClient(tc)
        fileSystem = new SeqeraFileSystem(this, datasetClient)
        fileSystem.registerHandler(new DatasetsResourceHandler(fileSystem, datasetClient))
        fileSystem.registerHandler(new DataLinksResourceHandler(fileSystem, new SeqeraDataLinkClient(tc)))
        return fileSystem
    }

    @Override
    synchronized FileSystem getFileSystem(URI uri) {
        checkScheme(uri)
        if (!fileSystem) throw new FileSystemNotFoundException("No seqera:// filesystem has been created yet")
        return fileSystem
    }

    synchronized SeqeraFileSystem getOrCreateFileSystem(URI uri, Map<String, ?> env) {
        checkScheme(uri)
        if (!fileSystem) newFileSystem(uri, env ?: Collections.<String, Object>emptyMap())
        return fileSystem
    }

    @Override
    SeqeraPath getPath(URI uri) {
        final fs = getOrCreateFileSystem(uri, Collections.emptyMap())
        return new SeqeraPath(fs, uri.toString())
    }

    // ---- read ----

    @Override
    InputStream newInputStream(Path path, OpenOption... options) throws IOException {
        final sp = toSeqeraPath(path)
        if (sp.depth() < 3)
            throw new IllegalArgumentException("newInputStream requires a leaf path: $path")
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final h = fs.getHandler(sp.resourceType)
        if (!h) throw new NoSuchFileException(path.toString(), null, "Unsupported resource type: ${sp.resourceType}")
        return h.newInputStream(sp)
    }

    @Override
    SeekableByteChannel newByteChannel(Path path, Set<? extends OpenOption> options, FileAttribute<?>... attrs) throws IOException {
        if (options?.contains(StandardOpenOption.WRITE) || options?.contains(StandardOpenOption.APPEND))
            throw new UnsupportedOperationException("seqera:// filesystem is read-only")
        return new DatasetInputStream(newInputStream(path))
    }

    // ---- attributes ----

    @Override
    <A extends BasicFileAttributes> A readAttributes(Path path, Class<A> type, LinkOption... options) throws IOException {
        if (!BasicFileAttributes.isAssignableFrom(type))
            throw new UnsupportedOperationException("Attribute type not supported: $type")
        final sp = toSeqeraPath(path)
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final d = sp.depth()
        if (d < 3) {
            validateSharedDirectoryExists(fs, sp)
            return (A) new SeqeraFileAttributes(true)
        }
        final h = fs.getHandler(sp.resourceType)
        if (!h) throw new NoSuchFileException(path.toString(), null, "Unsupported resource type: ${sp.resourceType}")
        return (A) h.readAttributes(sp)
    }

    @Override
    Map<String, Object> readAttributes(Path path, String attributes, LinkOption... options) throws IOException {
        throw new UnsupportedOperationException("readAttributes(String) not supported by `seqera://` filesystem")
    }

    // ---- access ----

    @Override
    void checkAccess(Path path, AccessMode... modes) throws IOException {
        final sp = toSeqeraPath(path)
        for (AccessMode m : modes) {
            if (m == AccessMode.WRITE || m == AccessMode.EXECUTE)
                throw new AccessDeniedException(path.toString(), null, "seqera:// filesystem is read-only")
        }
        final d = sp.depth()
        if (d == 0) return
        if (d < 3) {
            validateSharedDirectoryExists(sp.getFileSystem() as SeqeraFileSystem, sp)
            return
        }
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final h = fs.getHandler(sp.resourceType)
        if (!h) throw new NoSuchFileException(path.toString(), null, "Unsupported resource type: ${sp.resourceType}")
        h.checkAccess(sp, modes)
    }

    // ---- directory stream ----

    @Override
    DirectoryStream<Path> newDirectoryStream(Path dir, DirectoryStream.Filter<? super Path> filter) throws IOException {
        final sp = toSeqeraPath(dir)
        final fs = sp.getFileSystem() as SeqeraFileSystem
        final d = sp.depth()
        List<Path> entries
        if (d == 0) {
            fs.loadOrgWorkspaceCache()
            entries = fs.listOrgNames().collect { String org -> sp.resolve(org) as Path }
        } else if (d == 1) {
            fs.loadOrgWorkspaceCache()
            entries = fs.listWorkspaceNames(sp.org).collect { String ws -> sp.resolve(ws) as Path }
        } else if (d == 2) {
            fs.resolveWorkspaceId(sp.org, sp.workspace) // validates existence
            entries = fs.getResourceTypes().collect { String rt -> sp.resolve(rt) as Path }
        } else {
            final h = fs.getHandler(sp.resourceType)
            if (!h) throw new NoSuchFileException(dir.toString(), null, "Unsupported resource type: ${sp.resourceType}")
            entries = h.list(sp)
        }

        final filtered = filter ? entries.findAll { Path p ->
            try { filter.accept(p) }
            catch (IOException e) { throw new DirectoryIteratorException(e) }
        } : entries

        return new DirectoryStream<Path>() {
            @Override Iterator<Path> iterator() { filtered.iterator() }
            @Override void close() {}
        }
    }

    // ---- copy ----

    @Override
    void copy(Path source, Path target, CopyOption... options) throws IOException {
        toSeqeraPath(source)
        if (target instanceof SeqeraPath)
            throw new UnsupportedOperationException("seqera:// filesystem is read-only")
        try (final InputStream is = newInputStream(source)) {
            Files.copy(is, target, options)
        }
    }

    // ---- unsupported mutations ----

    @Override void move(Path s, Path t, CopyOption... o) { throw new UnsupportedOperationException("move() not supported") }
    @Override void delete(Path p) { throw new UnsupportedOperationException("delete() not supported") }
    @Override void createDirectory(Path d, FileAttribute<?>... a) { throw new UnsupportedOperationException("createDirectory() not supported") }
    @Override boolean isSameFile(Path a, Path b) { a == b }
    @Override boolean isHidden(Path p) { false }
    @Override FileStore getFileStore(Path p) { throw new UnsupportedOperationException("getFileStore() not supported") }
    @Override <V extends FileAttributeView> V getFileAttributeView(Path p, Class<V> t, LinkOption... o) { null }
    @Override void setAttribute(Path p, String a, Object v, LinkOption... o) { throw new UnsupportedOperationException("setAttribute() not supported") }

    // ---- helpers ----

    private static SeqeraPath toSeqeraPath(Path path) {
        if (path !instanceof SeqeraPath) throw new ProviderMismatchException()
        return (SeqeraPath) path
    }

    private static void checkScheme(URI uri) {
        if (uri.scheme?.toLowerCase() != SCHEME)
            throw new IllegalArgumentException("Not a seqera:// URI: $uri")
    }

    private static void validateSharedDirectoryExists(SeqeraFileSystem fs, SeqeraPath sp) throws NoSuchFileException {
        final d = sp.depth()
        if (d == 0) return
        fs.loadOrgWorkspaceCache()
        if (d >= 1 && !fs.listOrgNames().contains(sp.org))
            throw new NoSuchFileException("seqera://${sp.org}", null, "Organisation not found")
        if (d >= 2) fs.resolveWorkspaceId(sp.org, sp.workspace)
    }
}
```

- [ ] **Step 2: Compile**

Run: `./gradlew :plugins:nf-tower:compileGroovy`
Expected: **one** remaining compile failure — `import io.seqera.tower.plugin.datalink.SeqeraDataLinkClient` and `import io.seqera.tower.plugin.fs.handler.DataLinksResourceHandler` are not yet defined. These are created in T008–T010; continue without committing yet.

- [ ] **Step 3: Temporarily stub the missing imports so compile passes and existing tests can run**

Create a minimal stub at `plugins/nf-tower/src/main/io/seqera/tower/plugin/datalink/SeqeraDataLinkClient.groovy`:

```groovy
package io.seqera.tower.plugin.datalink

import groovy.transform.CompileStatic
import io.seqera.tower.plugin.TowerClient

@CompileStatic
class SeqeraDataLinkClient {
    SeqeraDataLinkClient(TowerClient tc) {}
}
```

Create a minimal stub at `plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/handler/DataLinksResourceHandler.groovy`:

```groovy
package io.seqera.tower.plugin.fs.handler

import java.nio.file.AccessMode
import java.nio.file.Path
import groovy.transform.CompileStatic
import io.seqera.tower.plugin.datalink.SeqeraDataLinkClient
import io.seqera.tower.plugin.fs.ResourceTypeHandler
import io.seqera.tower.plugin.fs.SeqeraFileAttributes
import io.seqera.tower.plugin.fs.SeqeraFileSystem
import io.seqera.tower.plugin.fs.SeqeraPath

@CompileStatic
class DataLinksResourceHandler implements ResourceTypeHandler {
    DataLinksResourceHandler(SeqeraFileSystem fs, SeqeraDataLinkClient client) {}
    @Override String getResourceType() { 'data-links' }
    @Override List<Path> list(SeqeraPath dir) { throw new UnsupportedOperationException('stub') }
    @Override SeqeraFileAttributes readAttributes(SeqeraPath path) { throw new UnsupportedOperationException('stub') }
    @Override InputStream newInputStream(SeqeraPath path) { throw new UnsupportedOperationException('stub') }
    @Override void checkAccess(SeqeraPath path, AccessMode... modes) { throw new UnsupportedOperationException('stub') }
}
```

- [ ] **Step 4: Compile and run the existing nf-tower tests**

Run: `./gradlew :plugins:nf-tower:test -i`
Expected: all existing tests pass — including `SeqeraDatasetClientTest`, `SeqeraFileSystemTest`, `SeqeraPathTest`, `SeqeraFileSystemProviderTest`, `DatasetsResourceHandlerTest` from T004. Any failure is a refactor regression and must be fixed before committing.

- [ ] **Step 5: Commit**

```bash
git add plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/SeqeraFileSystemProvider.groovy \
        plugins/nf-tower/src/main/io/seqera/tower/plugin/datalink/SeqeraDataLinkClient.groovy \
        plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/handler/DataLinksResourceHandler.groovy
git commit -s -m "refactor(nf-tower): SeqeraFileSystemProvider dispatches to ResourceTypeHandler; add data-link stubs"
```

**Checkpoint**: Refactor done. Dataset behavior is unchanged, routed through `DatasetsResourceHandler`. Data-link stubs exist but throw `UnsupportedOperationException`.

---

## Phase 2: Data-Link API Client

### T007 — Implement `SeqeraDataLinkClient`

**Files:**
- Modify: `plugins/nf-tower/src/main/io/seqera/tower/plugin/datalink/SeqeraDataLinkClient.groovy`

- [ ] **Step 1: Replace the stub with the real client**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * (Apache 2.0)
 */

package io.seqera.tower.plugin.datalink

import java.nio.file.AccessDeniedException
import java.nio.file.NoSuchFileException

import groovy.json.JsonSlurper
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.model.DataLinkContentResponse
import io.seqera.tower.model.DataLinkDownloadUrlResponse
import io.seqera.tower.model.DataLinkDto
import io.seqera.tower.model.DataLinkItem
import io.seqera.tower.model.DataLinkItemType
import io.seqera.tower.model.DataLinkProvider
import io.seqera.tower.plugin.TowerClient
import io.seqera.tower.plugin.exception.ForbiddenException
import io.seqera.tower.plugin.exception.NotFoundException
import io.seqera.tower.plugin.exception.UnauthorizedException
import nextflow.exception.AbortOperationException

/**
 * Typed client for Seqera Platform data-link API endpoints.
 * Delegates HTTP execution to {@link TowerClient#sendApiRequest}.
 */
@Slf4j
@CompileStatic
class SeqeraDataLinkClient {

    private static final int LIST_PAGE_SIZE = 100

    private final TowerClient towerClient

    SeqeraDataLinkClient(TowerClient tc) { this.towerClient = tc }

    private String getEndpoint() { towerClient.endpoint }

    /**
     * GET /data-links?workspaceId={ws}&max={n}&offset={o}
     * Exhausts pagination and returns all data-links in the workspace.
     */
    List<DataLinkDto> listDataLinks(long workspaceId) {
        final out = new ArrayList<DataLinkDto>()
        int offset = 0
        while (true) {
            final url = "${endpoint}/data-links?workspaceId=${workspaceId}&max=${LIST_PAGE_SIZE}&offset=${offset}"
            log.debug "SeqeraDataLinkClient GET $url"
            final resp = towerClient.sendApiRequest(url)
            checkFsResponse(resp, url)
            final json = new JsonSlurper().parseText(resp.message) as Map
            final items = json.dataLinks as List<Map>
            if (items) {
                for (Map m : items) out.add(mapDataLink(m))
                offset += items.size()
            }
            final total = (json.totalSize as Long) ?: 0L
            if (!items || offset >= total) break
        }
        return out
    }

    /**
     * GET /data-links/{id}/content?workspaceId={ws}&path={sub}&nextPageToken={tok}
     * Works for directories and single files. Exhausts {@code nextPageToken}.
     * Returns a synthesised {@link DataLinkContentResponse} with concatenated objects.
     */
    DataLinkContentResponse getContent(String dataLinkId, String subPath, long workspaceId) {
        final out = new DataLinkContentResponse()
        out.objects = new ArrayList<DataLinkItem>()
        String token = null
        while (true) {
            String url = "${endpoint}/data-links/${dataLinkId}/content?workspaceId=${workspaceId}"
            if (subPath) url += "&path=${encode(subPath)}"
            if (token) url += "&nextPageToken=${encode(token)}"
            log.debug "SeqeraDataLinkClient GET $url"
            final resp = towerClient.sendApiRequest(url)
            checkFsResponse(resp, url)
            final json = new JsonSlurper().parseText(resp.message) as Map
            if (out.originalPath == null) out.originalPath = json.originalPath as String
            final items = json.objects as List<Map>
            if (items) for (Map m : items) out.objects.add(mapItem(m))
            token = json.nextPageToken as String
            if (!token) break
        }
        return out
    }

    /** GET /data-links/{id}/download?workspaceId={ws}&path={sub} */
    DataLinkDownloadUrlResponse getDownloadUrl(String dataLinkId, String subPath, long workspaceId) {
        final url = "${endpoint}/data-links/${dataLinkId}/download?workspaceId=${workspaceId}&path=${encode(subPath ?: '')}"
        log.debug "SeqeraDataLinkClient GET $url"
        final resp = towerClient.sendApiRequest(url)
        checkFsResponse(resp, url)
        final json = new JsonSlurper().parseText(resp.message) as Map
        final out = new DataLinkDownloadUrlResponse()
        out.url = json.url as String
        return out
    }

    // ---- helpers ----

    private static String encode(String s) {
        new URI(null, null, s, null).rawPath
    }

    private static void checkFsResponse(TowerClient.Response resp, String url) {
        if (!resp.error) return
        final code = resp.code
        if (code == 401)
            throw new AbortOperationException("Seqera authentication failed — check tower.accessToken or TOWER_ACCESS_TOKEN")
        if (code == 403)
            throw new AccessDeniedException(url, null, "Forbidden — check workspace permissions")
        if (code == 404)
            throw new NoSuchFileException(url)
        throw new IOException("Seqera API error: HTTP ${code} for ${url}")
    }

    private static DataLinkDto mapDataLink(Map m) {
        final dto = new DataLinkDto()
        dto.id = m.id as String
        dto.name = m.name as String
        dto.description = m.description as String
        dto.resourceRef = m.resourceRef as String
        if (m.provider) dto.provider = DataLinkProvider.fromValue(m.provider as String)
        dto.region = m.region as String
        return dto
    }

    private static DataLinkItem mapItem(Map m) {
        final it = new DataLinkItem()
        it.name = m.name as String
        if (m.type) it.type = DataLinkItemType.fromValue(m.type as String)
        it.size = (m.size as Long) ?: 0L
        it.mimeType = m.mimeType as String
        return it
    }
}
```

- [ ] **Step 2: Verify `DataLinkProvider.fromValue` and `DataLinkItemType.fromValue` exist**

Run: `javap -p /home/jorgee/IdeaProjects/nextflow/plugins/nf-tower/build/target/libs/tower-api-1.121.0.jar | grep -A2 'class io.seqera.tower.model.DataLinkProvider' | head -20` — or just proceed. These `fromValue` methods are standard on generated Micronaut enums; if not present the compile step will tell us.

- [ ] **Step 3: Compile**

Run: `./gradlew :plugins:nf-tower:compileGroovy`
Expected: success. If `fromValue` is missing, fall back to `DataLinkProvider.values().find { it.toString() == m.provider }`.

- [ ] **Step 4: Commit**

```bash
git add plugins/nf-tower/src/main/io/seqera/tower/plugin/datalink/SeqeraDataLinkClient.groovy
git commit -s -m "feat(nf-tower): add SeqeraDataLinkClient with pagination and error mapping"
```

### T008 — Unit tests for `SeqeraDataLinkClient`

**Files:**
- Create: `plugins/nf-tower/src/test/io/seqera/tower/plugin/datalink/SeqeraDataLinkClientTest.groovy`

- [ ] **Step 1: Write the Spock spec**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * (Apache 2.0)
 */
package io.seqera.tower.plugin.datalink

import java.nio.file.AccessDeniedException
import java.nio.file.NoSuchFileException

import groovy.json.JsonOutput
import io.seqera.tower.plugin.TowerClient
import nextflow.exception.AbortOperationException
import spock.lang.Specification

class SeqeraDataLinkClientTest extends Specification {

    private static final String EP = 'https://api.example.com'

    private TowerClient tower() {
        def tc = Spy(TowerClient)
        tc.@endpoint = EP
        return tc
    }

    private static TowerClient.Response ok(String body) { new TowerClient.Response(200, body) }
    private static TowerClient.Response err(int code)    { new TowerClient.Response(code, "error $code") }

    // ---- listDataLinks ----

    def "listDataLinks returns parsed DTOs for a single page"() {
        given:
        def body = JsonOutput.toJson([dataLinks: [
            [id: 'dl-1', name: 'inputs', provider: 'aws', resourceRef: 's3://bucket/'],
            [id: 'dl-2', name: 'archive', provider: 'google', resourceRef: 'gs://bucket/']
        ], totalSize: 2])
        def tc = tower()
        tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=0") >> ok(body)
        def client = new SeqeraDataLinkClient(tc)

        when:
        def list = client.listDataLinks(10L)

        then:
        list.size() == 2
        list[0].id == 'dl-1'
        list[0].name == 'inputs'
        list[1].provider.toString() == 'google'
    }

    def "listDataLinks exhausts pagination"() {
        given:
        def page1 = JsonOutput.toJson([dataLinks: [[id: 'dl-1', name: 'a', provider: 'aws']], totalSize: 3])
        def page2 = JsonOutput.toJson([dataLinks: [[id: 'dl-2', name: 'b', provider: 'aws']], totalSize: 3])
        def page3 = JsonOutput.toJson([dataLinks: [[id: 'dl-3', name: 'c', provider: 'aws']], totalSize: 3])
        def tc = tower()
        tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=0") >> ok(page1)
        tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=1") >> ok(page2)
        tc.sendApiRequest("${EP}/data-links?workspaceId=10&max=100&offset=2") >> ok(page3)
        def client = new SeqeraDataLinkClient(tc)

        when:
        def list = client.listDataLinks(10L)

        then:
        list*.id == ['dl-1', 'dl-2', 'dl-3']
    }

    // ---- getContent ----

    def "getContent single page returns parsed items"() {
        given:
        def body = JsonOutput.toJson([
            originalPath: 'reads/',
            objects: [
                [name: 'a.fq', type: 'FILE', size: 123, mimeType: 'application/gzip'],
                [name: 'b.fq', type: 'FILE', size: 456, mimeType: 'application/gzip']
            ]])
        def tc = tower()
        tc.sendApiRequest("${EP}/data-links/dl-1/content?workspaceId=10&path=reads/") >> ok(body)
        def client = new SeqeraDataLinkClient(tc)

        when:
        def resp = client.getContent('dl-1', 'reads/', 10L)

        then:
        resp.objects.size() == 2
        resp.objects[0].name == 'a.fq'
        resp.objects[0].size == 123L
        resp.objects[0].type.toString() == 'FILE'
    }

    def "getContent follows nextPageToken"() {
        given:
        def p1 = JsonOutput.toJson([originalPath: '', objects: [[name: 'a', type: 'FILE', size: 1]], nextPageToken: 'T2'])
        def p2 = JsonOutput.toJson([originalPath: '', objects: [[name: 'b', type: 'FILE', size: 2]]])
        def tc = tower()
        tc.sendApiRequest("${EP}/data-links/dl-1/content?workspaceId=10") >> ok(p1)
        tc.sendApiRequest("${EP}/data-links/dl-1/content?workspaceId=10&nextPageToken=T2") >> ok(p2)
        def client = new SeqeraDataLinkClient(tc)

        when:
        def resp = client.getContent('dl-1', null, 10L)

        then:
        resp.objects*.name == ['a', 'b']
    }

    // ---- getDownloadUrl ----

    def "getDownloadUrl returns the signed URL"() {
        given:
        def tc = tower()
        tc.sendApiRequest("${EP}/data-links/dl-1/download?workspaceId=10&path=reads/a.fq") >> ok(JsonOutput.toJson([url: 'https://signed']))
        def client = new SeqeraDataLinkClient(tc)

        when:
        def resp = client.getDownloadUrl('dl-1', 'reads/a.fq', 10L)

        then:
        resp.url == 'https://signed'
    }

    // ---- error mapping ----

    def "401 throws AbortOperationException"() {
        given:
        def tc = tower()
        tc.sendApiRequest(_) >> err(401)
        def client = new SeqeraDataLinkClient(tc)

        when:
        client.listDataLinks(10L)

        then:
        thrown(AbortOperationException)
    }

    def "403 throws AccessDeniedException"() {
        given:
        def tc = tower()
        tc.sendApiRequest(_) >> err(403)
        def client = new SeqeraDataLinkClient(tc)

        when:
        client.getContent('dl-1', '', 10L)

        then:
        thrown(AccessDeniedException)
    }

    def "404 throws NoSuchFileException"() {
        given:
        def tc = tower()
        tc.sendApiRequest(_) >> err(404)
        def client = new SeqeraDataLinkClient(tc)

        when:
        client.getDownloadUrl('dl-1', 'missing', 10L)

        then:
        thrown(NoSuchFileException)
    }
}
```

- [ ] **Step 2: Run the tests**

Run: `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.datalink.SeqeraDataLinkClientTest' -i`
Expected: all pass.

- [ ] **Step 3: Commit**

```bash
git add plugins/nf-tower/src/test/io/seqera/tower/plugin/datalink/SeqeraDataLinkClientTest.groovy
git commit -s -m "test(nf-tower): SeqeraDataLinkClient unit tests — pagination, parse, error mapping"
```

---

## Phase 3: User Story 1 — Read a File Inside a Data-Link (P1) 🎯 MVP

**Goal**: `file('seqera://<org>/<ws>/data-links/<prov>/<name>/<path>')` returns an `InputStream` over the file content.

### T009 [US1] — Implement `DataLinksResourceHandler` (real, replacing the stub)

**Files:**
- Modify: `plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/handler/DataLinksResourceHandler.groovy`

- [ ] **Step 1: Replace the stub with the full handler**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * (Apache 2.0)
 */

package io.seqera.tower.plugin.fs.handler

import java.net.http.HttpClient
import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.file.AccessDeniedException
import java.nio.file.AccessMode
import java.nio.file.NoSuchFileException
import java.nio.file.Path
import java.time.Duration
import java.time.Instant

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.tower.model.DataLinkContentResponse
import io.seqera.tower.model.DataLinkDto
import io.seqera.tower.model.DataLinkItem
import io.seqera.tower.model.DataLinkItemType
import io.seqera.tower.plugin.datalink.SeqeraDataLinkClient
import io.seqera.tower.plugin.fs.ResourceTypeHandler
import io.seqera.tower.plugin.fs.SeqeraFileAttributes
import io.seqera.tower.plugin.fs.SeqeraFileSystem
import io.seqera.tower.plugin.fs.SeqeraPath

/**
 * {@link ResourceTypeHandler} for {@code data-links} resource type.
 * Listings and attributes go through the Seqera Platform API; file reads go through
 * a pre-signed URL fetched with a plain JDK {@link HttpClient} — the signed URL is
 * addressed to the cloud backend and the Seqera {@code Authorization} header must not be sent.
 */
@Slf4j
@CompileStatic
class DataLinksResourceHandler implements ResourceTypeHandler {

    public static final String TYPE = 'data-links'

    private final SeqeraFileSystem fs
    private final SeqeraDataLinkClient client
    private final HttpClient httpClient

    /** workspaceId → data-link list */
    private final Map<Long, List<DataLinkDto>> dataLinkCache = new LinkedHashMap<>()

    DataLinksResourceHandler(SeqeraFileSystem fs, SeqeraDataLinkClient client) {
        this(fs, client, HttpClient.newBuilder().connectTimeout(Duration.ofSeconds(10)).build())
    }

    /** Test-only constructor to inject a mock {@link HttpClient}. */
    DataLinksResourceHandler(SeqeraFileSystem fs, SeqeraDataLinkClient client, HttpClient httpClient) {
        this.fs = fs
        this.client = client
        this.httpClient = httpClient
    }

    @Override String getResourceType() { TYPE }

    @Override
    List<Path> list(SeqeraPath dir) throws IOException {
        final workspaceId = fs.resolveWorkspaceId(dir.org, dir.workspace)
        final trail = dir.trail
        if (trail.isEmpty()) {
            // List distinct providers in use
            final providers = resolveDataLinks(workspaceId)
                    .collect { DataLinkDto dl -> dl.provider?.toString() }
                    .findAll { String p -> p }
                    .toSet()
            return providers.toList().sort().collect { String p -> dir.resolve(p) as Path }
        }
        if (trail.size() == 1) {
            // List data-link names under the given provider
            final prov = trail[0]
            final names = resolveDataLinks(workspaceId)
                    .findAll { DataLinkDto dl -> dl.provider?.toString() == prov }
                    .collect { DataLinkDto dl -> dl.name }
                    .sort()
            if (names.isEmpty())
                throw new NoSuchFileException(dir.toString(), null, "No data-links for provider '$prov' in workspace '${dir.workspace}'")
            return names.collect { String n -> dir.resolve(n) as Path }
        }
        // trail.size() >= 2 — browse inside the data-link
        final dl = resolveDataLink(workspaceId, trail[0], trail[1])
        final subPath = trail.size() > 2 ? trail.subList(2, trail.size()).join('/') : ''
        final resp = client.getContent(dl.id, subPath, workspaceId)
        return (resp.objects ?: []).collect { DataLinkItem it -> dir.resolve(it.name) as Path }
    }

    @Override
    SeqeraFileAttributes readAttributes(SeqeraPath p) throws IOException {
        final workspaceId = fs.resolveWorkspaceId(p.org, p.workspace)
        final trail = p.trail
        // data-links/ dir itself, provider dir, data-link root — all directories
        if (trail.size() < 2) return new SeqeraFileAttributes(true)
        final dl = resolveDataLink(workspaceId, trail[0], trail[1])
        if (trail.size() == 2) return new SeqeraFileAttributes(true) // data-link root
        final subPath = trail.subList(2, trail.size()).join('/')
        final resp = client.getContent(dl.id, subPath, workspaceId)
        return attributesFor(resp, subPath, p)
    }

    @Override
    InputStream newInputStream(SeqeraPath p) throws IOException {
        if (p.trail.size() < 3)
            throw new IllegalArgumentException("newInputStream requires a file path inside a data-link: $p")
        final workspaceId = fs.resolveWorkspaceId(p.org, p.workspace)
        final dl = resolveDataLink(workspaceId, p.trail[0], p.trail[1])
        final subPath = p.trail.subList(2, p.trail.size()).join('/')
        final urlResp = client.getDownloadUrl(dl.id, subPath, workspaceId)
        if (!urlResp.url)
            throw new NoSuchFileException(p.toString(), null, "Platform returned no download URL")
        return fetchSignedUrl(urlResp.url)
    }

    @Override
    void checkAccess(SeqeraPath p, AccessMode... modes) throws IOException {
        for (AccessMode m : modes) {
            if (m == AccessMode.WRITE || m == AccessMode.EXECUTE)
                throw new AccessDeniedException(p.toString(), null, "seqera:// data-links are read-only")
        }
        // READ: rely on readAttributes to validate existence
        readAttributes(p)
    }

    // ---- private ----

    private synchronized List<DataLinkDto> resolveDataLinks(long workspaceId) {
        def cached = dataLinkCache.get(workspaceId)
        if (cached == null) {
            cached = client.listDataLinks(workspaceId)
            dataLinkCache.put(workspaceId, cached)
        }
        return cached
    }

    private DataLinkDto resolveDataLink(long workspaceId, String provider, String name) throws NoSuchFileException {
        final found = resolveDataLinks(workspaceId).find { DataLinkDto dl ->
            dl.provider?.toString() == provider && dl.name == name
        }
        if (!found)
            throw new NoSuchFileException(
                "seqera://.../data-links/${provider}/${name}",
                null,
                "Data-link '${name}' not found for provider '${provider}' in workspace '$workspaceId'")
        return found
    }

    private SeqeraFileAttributes attributesFor(DataLinkContentResponse resp, String subPath, SeqeraPath pathForErrors) throws NoSuchFileException {
        final items = resp.objects ?: []
        // Single-file content response: one object whose name matches the final segment
        final lastSeg = subPath.contains('/') ? subPath.substring(subPath.lastIndexOf('/') + 1) : subPath
        final single = items.find { DataLinkItem it -> it.name == lastSeg && it.type == DataLinkItemType.FILE }
        if (single)
            return new SeqeraFileAttributes(single.size ?: 0L, Instant.EPOCH, Instant.EPOCH, pathForErrors.toString())
        // Otherwise treat as a directory (content response with multiple children or zero)
        // If there are no children AND no originalPath, the path does not exist
        if (items.isEmpty() && !resp.originalPath)
            throw new NoSuchFileException(pathForErrors.toString(), null, "Path not found inside data-link")
        return new SeqeraFileAttributes(true)
    }

    private InputStream fetchSignedUrl(String url) throws IOException {
        final req = HttpRequest.newBuilder(URI.create(url))
                .timeout(Duration.ofMinutes(5))
                .GET()
                .build()
        try {
            final HttpResponse<InputStream> resp = httpClient.send(req, HttpResponse.BodyHandlers.ofInputStream())
            final status = resp.statusCode()
            if (status >= 200 && status < 300) return resp.body()
            resp.body()?.close()
            throw new IOException("Signed URL fetch failed: HTTP $status for $url")
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt()
            throw new IOException("Interrupted while fetching signed URL", e)
        }
    }
}
```

- [ ] **Step 2: Compile**

Run: `./gradlew :plugins:nf-tower:compileGroovy`
Expected: success.

- [ ] **Step 3: Commit**

```bash
git add plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/handler/DataLinksResourceHandler.groovy
git commit -s -m "feat(nf-tower): implement DataLinksResourceHandler (list, readAttributes, newInputStream)"
```

### T010 [US1] — Unit tests for `DataLinksResourceHandler.newInputStream`

**Files:**
- Create: `plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/handler/DataLinksResourceHandlerTest.groovy`

- [ ] **Step 1: Write the spec — MVP scenarios for newInputStream**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * (Apache 2.0)
 */
package io.seqera.tower.plugin.fs.handler

import java.net.http.HttpClient
import java.net.http.HttpResponse
import java.nio.file.NoSuchFileException

import io.seqera.tower.model.DataLinkContentResponse
import io.seqera.tower.model.DataLinkDownloadUrlResponse
import io.seqera.tower.model.DataLinkDto
import io.seqera.tower.model.DataLinkItem
import io.seqera.tower.model.DataLinkItemType
import io.seqera.tower.model.DataLinkProvider
import io.seqera.tower.plugin.datalink.SeqeraDataLinkClient
import io.seqera.tower.plugin.fs.SeqeraFileSystem
import io.seqera.tower.plugin.fs.SeqeraPath
import spock.lang.Specification

class DataLinksResourceHandlerTest extends Specification {

    private SeqeraFileSystem fs = Mock(SeqeraFileSystem)
    private SeqeraDataLinkClient client = Mock(SeqeraDataLinkClient)
    private HttpClient http = Mock(HttpClient)
    private DataLinksResourceHandler handler = new DataLinksResourceHandler(fs, client, http)

    private DataLinkDto dl(String id, String name, DataLinkProvider p) {
        def d = new DataLinkDto(); d.id = id; d.name = name; d.provider = p; return d
    }
    private DataLinkItem item(String name, DataLinkItemType t, long size) {
        def i = new DataLinkItem(); i.name = name; i.type = t; i.size = size; return i
    }

    // ---- newInputStream ----

    def "newInputStream resolves (provider,name,subPath) and streams the signed URL"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/reads/a.fq')
        def signedBody = new ByteArrayInputStream('data'.bytes)
        def httpResp = Mock(HttpResponse) {
            statusCode() >> 200
            body() >> signedBody
        }
        def urlResp = new DataLinkDownloadUrlResponse(); urlResp.url = 'https://signed/a'

        when:
        def stream = handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'inputs', DataLinkProvider.AWS)]
        1 * client.getDownloadUrl('dl-1', 'reads/a.fq', 10L) >> urlResp
        1 * http.send(_, _) >> httpResp
        stream === signedBody
    }

    def "newInputStream throws NoSuchFileException when data-link unknown"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/unknown/reads/a.fq')

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'inputs', DataLinkProvider.AWS)]
        thrown(NoSuchFileException)
    }

    def "newInputStream requires trail.size >= 3 (file path, not the data-link root itself)"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs')

        when:
        handler.newInputStream(path)

        then:
        thrown(IllegalArgumentException)
    }

    def "newInputStream wraps signed-URL HTTP 403 as IOException"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/reads/a.fq')
        def urlResp = new DataLinkDownloadUrlResponse(); urlResp.url = 'https://signed/a'
        def httpResp = Mock(HttpResponse) {
            statusCode() >> 403
            body() >> new ByteArrayInputStream(new byte[0])
        }

        when:
        handler.newInputStream(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'inputs', DataLinkProvider.AWS)]
        1 * client.getDownloadUrl('dl-1', 'reads/a.fq', 10L) >> urlResp
        1 * http.send(_, _) >> httpResp
        thrown(IOException)
    }
}
```

- [ ] **Step 2: Run the tests**

Run: `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.handler.DataLinksResourceHandlerTest' -i`
Expected: all pass.

- [ ] **Step 3: Commit**

```bash
git add plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/handler/DataLinksResourceHandlerTest.groovy
git commit -s -m "test(nf-tower): DataLinksResourceHandler.newInputStream unit tests"
```

**Checkpoint**: US1 complete — file reads through data-link paths work end to end in unit tests.

---

## Phase 4: User Story 2 — Browse Data-Link Hierarchy (P2)

### T011 [US2] — List & readAttributes tests for DataLinksResourceHandler

**Files:**
- Modify: `plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/handler/DataLinksResourceHandlerTest.groovy`

- [ ] **Step 1: Append list / readAttributes tests**

Add the following specs to `DataLinksResourceHandlerTest`:

```groovy
    // ---- list: depth 3 (data-links/) ----

    def "list at data-links/ returns distinct providers in use"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links')

        when:
        def paths = handler.list(path)

        then:
        1 * fs.resolveWorkspaceId('acme', 'research') >> 10L
        1 * client.listDataLinks(10L) >> [
                dl('dl-1', 'a', DataLinkProvider.AWS),
                dl('dl-2', 'b', DataLinkProvider.GOOGLE),
                dl('dl-3', 'c', DataLinkProvider.AWS)
        ]
        paths*.toString().sort() == [
                'seqera://acme/research/data-links/AWS',
                'seqera://acme/research/data-links/GOOGLE'
        ]
    }

    // ---- list: depth 4 (data-links/<prov>/) ----

    def "list at data-links/<provider>/ returns data-link names for that provider"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS')

        when:
        def paths = handler.list(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> [
                dl('dl-1', 'inputs', DataLinkProvider.AWS),
                dl('dl-2', 'archive', DataLinkProvider.AWS),
                dl('dl-3', 'onGcs', DataLinkProvider.GOOGLE)
        ]
        paths*.toString() == [
                'seqera://acme/research/data-links/AWS/archive',
                'seqera://acme/research/data-links/AWS/inputs'
        ]
    }

    def "list at data-links/<provider>/ throws when no data-links for that provider"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AZURE')

        when:
        handler.list(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'x', DataLinkProvider.AWS)]
        thrown(NoSuchFileException)
    }

    // ---- list: depth 5 (data-link root) ----

    def "list at data-link root returns top-level objects"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs')
        def content = new DataLinkContentResponse()
        content.objects = [item('reads', DataLinkItemType.FOLDER, 0), item('samplesheet.csv', DataLinkItemType.FILE, 42)]

        when:
        def paths = handler.list(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'inputs', DataLinkProvider.AWS)]
        1 * client.getContent('dl-1', '', 10L) >> content
        paths*.toString() == [
                'seqera://acme/research/data-links/AWS/inputs/reads',
                'seqera://acme/research/data-links/AWS/inputs/samplesheet.csv'
        ]
    }

    // ---- list: depth 6+ (nested sub-path) ----

    def "list at deep sub-path browses the correct sub-path"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/reads')
        def content = new DataLinkContentResponse()
        content.objects = [item('a.fq', DataLinkItemType.FILE, 1), item('b.fq', DataLinkItemType.FILE, 2)]

        when:
        def paths = handler.list(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'inputs', DataLinkProvider.AWS)]
        1 * client.getContent('dl-1', 'reads', 10L) >> content
        paths*.toString() == [
                'seqera://acme/research/data-links/AWS/inputs/reads/a.fq',
                'seqera://acme/research/data-links/AWS/inputs/reads/b.fq'
        ]
    }

    // ---- readAttributes ----

    def "readAttributes at data-links/ resource-type dir reports directory"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        attr.directory
        !attr.regularFile
    }

    def "readAttributes at data-link root reports directory"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs')

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'inputs', DataLinkProvider.AWS)]
        attr.directory
    }

    def "readAttributes on a file sub-path reports file with size"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/reads/a.fq')
        def content = new DataLinkContentResponse()
        content.objects = [item('a.fq', DataLinkItemType.FILE, 123)]

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'inputs', DataLinkProvider.AWS)]
        1 * client.getContent('dl-1', 'reads/a.fq', 10L) >> content
        attr.regularFile
        attr.size() == 123L
    }

    def "readAttributes on a directory sub-path reports directory"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/reads')
        def content = new DataLinkContentResponse()
        content.originalPath = 'reads/'
        content.objects = [item('a.fq', DataLinkItemType.FILE, 1), item('b.fq', DataLinkItemType.FILE, 2)]

        when:
        def attr = handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'inputs', DataLinkProvider.AWS)]
        1 * client.getContent('dl-1', 'reads', 10L) >> content
        attr.directory
    }
```

- [ ] **Step 2: Run the tests**

Run: `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.handler.DataLinksResourceHandlerTest' -i`
Expected: all pass.

- [ ] **Step 3: Commit**

```bash
git add plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/handler/DataLinksResourceHandlerTest.groovy
git commit -s -m "test(nf-tower): DataLinksResourceHandler list & readAttributes unit tests"
```

### T012 [US2] — Provider-level browsing test (workspace listing enumerates handlers)

**Files:**
- Modify: `plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/SeqeraFileSystemProviderTest.groovy`

- [ ] **Step 1: Add a spec confirming `datasets` AND `data-links` appear when listing a workspace**

Append to the test class (before the final `}`):

```groovy
    def "listing a workspace enumerates both registered resource types"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())

        final provider = new SeqeraFileSystemProvider()
        provider.newFileSystem(URI.create('seqera://'), [:])
        final fs = provider.getFileSystem(URI.create('seqera://')) as SeqeraFileSystem
        final wsPath = new SeqeraPath(fs, 'seqera://acme/research')

        when:
        final List<java.nio.file.Path> entries = []
        provider.newDirectoryStream(wsPath, null).withCloseable { s -> s.iterator().each { entries.add(it) } }

        then:
        entries*.toString().sort() == [
                'seqera://acme/research/data-links',
                'seqera://acme/research/datasets'
        ]
    }
```

- [ ] **Step 2: Run the test**

Run: `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.SeqeraFileSystemProviderTest' -i`
Expected: all pass — includes both the new spec and the existing dataset-related specs.

- [ ] **Step 3: Commit**

```bash
git add plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/SeqeraFileSystemProviderTest.groovy
git commit -s -m "test(nf-tower): workspace listing enumerates both datasets and data-links"
```

**Checkpoint**: US2 complete — browsing works at every depth.

---

## Phase 5: User Story 3 — Meaningful Errors (P3)

### T013 [US3] — Error-mapping tests for data-link paths

**Files:**
- Modify: `plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/handler/DataLinksResourceHandlerTest.groovy`

- [ ] **Step 1: Append error-path specs**

```groovy
    def "unknown data-link under a known provider throws with clear message"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/ghost/a.fq')

        when:
        handler.list(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'inputs', DataLinkProvider.AWS)]
        def ex = thrown(NoSuchFileException)
        ex.message.toLowerCase().contains('not found') || ex.reason?.toLowerCase()?.contains('not found')
    }

    def "missing sub-path inside a data-link surfaces as NoSuchFileException"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/does/not/exist')
        def empty = new DataLinkContentResponse()
        empty.originalPath = null
        empty.objects = []

        when:
        handler.readAttributes(path)

        then:
        1 * fs.resolveWorkspaceId(_, _) >> 10L
        1 * client.listDataLinks(10L) >> [dl('dl-1', 'inputs', DataLinkProvider.AWS)]
        1 * client.getContent('dl-1', 'does/not/exist', 10L) >> empty
        thrown(NoSuchFileException)
    }

    def "checkAccess with WRITE is rejected"() {
        given:
        def path = new SeqeraPath(fs, 'seqera://acme/research/data-links/AWS/inputs/a.fq')

        when:
        handler.checkAccess(path, java.nio.file.AccessMode.WRITE)

        then:
        thrown(java.nio.file.AccessDeniedException)
    }
```

- [ ] **Step 2: Run the tests**

Run: `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.handler.DataLinksResourceHandlerTest' -i`
Expected: all pass.

- [ ] **Step 3: Commit**

```bash
git add plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/handler/DataLinksResourceHandlerTest.groovy
git commit -s -m "test(nf-tower): data-link error paths — unknown link, missing sub-path, write-rejected"
```

### T014 [US3] — Unsupported-resource-type error via the provider

**Files:**
- Modify: `plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/SeqeraFileSystemProviderTest.groovy`

- [ ] **Step 1: Add a provider-level dispatch test**

Append:

```groovy
    def "newInputStream on an unsupported resource type throws NoSuchFileException"() {
        given:
        def tc = spyTower()
        tc.sendApiRequest("${ENDPOINT}/user-info") >> ok(userInfoJson())
        tc.sendApiRequest("${ENDPOINT}/user/42/workspaces") >> ok(workspacesJson())
        def provider = new SeqeraFileSystemProvider()
        provider.newFileSystem(URI.create('seqera://'), [:])
        def fs = provider.getFileSystem(URI.create('seqera://')) as SeqeraFileSystem
        def path = new SeqeraPath(fs, 'seqera://acme/research/unknown-type/foo')

        when:
        provider.newInputStream(path)

        then:
        def ex = thrown(NoSuchFileException)
        ex.message.contains('Unsupported resource type') || ex.reason?.contains('Unsupported resource type')
    }
```

- [ ] **Step 2: Run the tests**

Run: `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.SeqeraFileSystemProviderTest' -i`
Expected: all pass.

- [ ] **Step 3: Commit**

```bash
git add plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/SeqeraFileSystemProviderTest.groovy
git commit -s -m "test(nf-tower): unsupported resource type surfaces NoSuchFileException"
```

**Checkpoint**: US3 complete — all error paths produce clear, type-specific exceptions.

---

## Phase 6: User Story 4 — Extensibility Validation (P4)

### T015 [US4] — Architectural guard test

**Files:**
- Create: `plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/ResourceTypeAbstractionTest.groovy`

- [ ] **Step 1: Write a guard test**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 * (Apache 2.0)
 */
package io.seqera.tower.plugin.fs

import spock.lang.Specification

/**
 * Guards that the generic NIO layer does not reach into resource-type-specific packages.
 * {@link SeqeraPath}, {@link SeqeraFileSystem}, {@link SeqeraFileSystemProvider} must not
 * depend on {@code dataset/}, {@code datalink/}, or {@code fs/handler/} — handlers dispatch,
 * but dispatch lives behind the {@link ResourceTypeHandler} interface.
 */
class ResourceTypeAbstractionTest extends Specification {

    static final Class[] GENERIC_CLASSES = [SeqeraPath, SeqeraFileSystem, SeqeraFileAttributes]

    def "generic fs classes do not import resource-type-specific packages"() {
        expect:
        GENERIC_CLASSES.each { Class c ->
            final src = new File("plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/${c.simpleName}.groovy").text
            assert !src.contains('io.seqera.tower.plugin.dataset.'), "${c.simpleName} must not import dataset package"
            assert !src.contains('io.seqera.tower.plugin.datalink.'), "${c.simpleName} must not import datalink package"
            assert !src.contains('io.seqera.tower.plugin.fs.handler.'), "${c.simpleName} must not import handler package"
            assert !src.contains('DataLink') , "${c.simpleName} must not reference data-link types"
            assert !src.contains('DatasetDto'), "${c.simpleName} must not reference DatasetDto"
        }
    }

    def "both handlers implement the ResourceTypeHandler interface"() {
        expect:
        ResourceTypeHandler.isAssignableFrom(io.seqera.tower.plugin.fs.handler.DatasetsResourceHandler)
        ResourceTypeHandler.isAssignableFrom(io.seqera.tower.plugin.fs.handler.DataLinksResourceHandler)
    }
}
```

- [ ] **Step 2: Run**

Run: `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.ResourceTypeAbstractionTest' -i`
Expected: all pass. If any import lingered from the refactor (e.g. `SeqeraFileSystem.groovy` still references `DatasetDto` via a cache field type), fix the import before proceeding — the dataset-specific caches belong in `DatasetsResourceHandler` long-term, but keeping the `DatasetDto` typed cache in `SeqeraFileSystem` **is acceptable** as long as the import is `io.seqera.tower.model.DatasetDto` (the tower-api DTO, not the plugin's `dataset` package). The test check against `'DatasetDto'` in source text guards the generic classes; if the existing field trips this, refactor the cache into the handler (move `datasetCache`, `versionCache`, `resolveDatasets`, `resolveDataset`, `resolveVersions`, `invalidateDatasetCache` from `SeqeraFileSystem` into `DatasetsResourceHandler`).

**Important refactor note**: if the guard test fails on `SeqeraFileSystem.groovy`, perform this sub-step:

- Remove `datasetCache`, `versionCache`, `resolveDatasets`, `resolveDataset`, `resolveVersions`, and `invalidateDatasetCache` from `SeqeraFileSystem.groovy`, along with their imports (`io.seqera.tower.model.DatasetDto`, `io.seqera.tower.model.DatasetVersionDto`).
- Move the same caches and methods into `DatasetsResourceHandler.groovy` as private fields and synchronized methods.
- Update `DatasetsResourceHandler` to use its own cache methods instead of calling `fs.resolveDatasets(...)` etc.
- Re-run the guard test to confirm.

- [ ] **Step 3: Commit**

```bash
git add plugins/nf-tower/src/test/io/seqera/tower/plugin/fs/ResourceTypeAbstractionTest.groovy \
        plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/SeqeraFileSystem.groovy \
        plugins/nf-tower/src/main/io/seqera/tower/plugin/fs/handler/DatasetsResourceHandler.groovy
git commit -s -m "test(nf-tower): enforce resource-type-agnostic boundary in generic fs classes"
```

**Checkpoint**: US4 complete — abstraction is validated by automated guard.

---

## Phase 7: Final Verification

**Note**: both the `plugins/nf-tower/VERSION` bump and `changelog.txt` entries are handled at release time by the repo's release process (see `CLAUDE.md § Release process`), not per-feature. This phase only verifies the build is green.

### T016 — Final full test run

- [ ] **Step 1: Run the full nf-tower test suite**

Run: `./gradlew :plugins:nf-tower:test -i`
Expected: all tests pass — dataset, data-link, path, provider, filesystem, abstraction guard.

- [ ] **Step 2: Run the full plugin build**

Run: `./gradlew :plugins:nf-tower:build`
Expected: success.

- [ ] **Step 3: Confirm no cloud-SDK dependencies were introduced**

Run:
```
./gradlew :plugins:nf-tower:dependencies --configuration runtimeClasspath | grep -iE 'aws-sdk|google-cloud|azure-' || echo 'OK: no cloud SDKs on classpath'
```
Expected output ends with `OK: no cloud SDKs on classpath`. (SC-006)

- [ ] **Step 4: Summary — nothing to commit; just confirm build-green at HEAD.**

---

## Appendix — Task Dependency Graph

```
T001 ─┐
T002 ─┼──► T003 ──► T004 ──► T005 ──► T006 ──► T007 ──► T008
      │                                              │
      │                                              └──► T009 ──► T010
      │                                                                │
      │                                                                ├──► T011
      │                                                                ├──► T012
      │                                                                ├──► T013
      │                                                                └──► T014
      │                                                                         │
      │                                                                         └──► T015
      │                                                                                  │
      └──────────────────────────────────────────────────────────────────────────────────┴──► T016
```
