# Implementation Plan: Seqera NIO Filesystem Support for Platform Data-Links

**Branch**: `260422-seqera-datalinks-fs` | **Date**: 2026-04-23 | **Spec**: [spec.md](spec.md) | **ADR**: [20260422-seqera-datalinks-filesystem](../../adr/20260422-seqera-datalinks-filesystem.md)

---

## Summary

Extend the `nf-tower` plugin's `seqera://` NIO filesystem (shipped in [260310-seqera-dataset-fs](../260310-seqera-dataset-fs/spec.md)) with a second resource type, `data-links`. Paths of the form `seqera://<org>/<ws>/data-links/<provider>/<name>/<sub-path>` resolve to entries inside Platform-managed data-links (S3/GCS/Azure buckets or prefixes). Listings and attribute queries hit the Platform's `/data-links/{id}/browse[/path]` endpoints; byte reads go through a pre-signed URL returned by `/data-links/{id}/generate-download-url` and fetched with a plain JDK HTTP client — no cloud SDK dependency.

As part of this change, extract a real `ResourceTypeHandler` abstraction from the existing dataset logic. `DatasetsResourceHandler` and `DataLinksResourceHandler` are parallel implementations; the generic classes (`SeqeraFileSystemProvider`, `SeqeraFileSystem`, `SeqeraPath`, `SeqeraFileAttributes`) become resource-type-agnostic for depth ≥ 3.

---

## Technical Context

**Language/Version**: Groovy 4.0.29, targeting Java 17 runtime
**Primary Dependencies**: `io.seqera:tower-api:1.121.0` (existing), `io.seqera:lib-httpx:2.1.0` (existing via TowerClient)
**Storage**: None (in-memory org/workspace + data-link list caches on `SeqeraFileSystem`)
**Testing**: Spock Framework with `Mock(TowerClient)` + `JsonOutput` fixtures — matches the dataset test style
**Target Platform**: nf-tower plugin; runs wherever Nextflow runs
**Performance Goals**: Listing any hierarchy level in ≤ 5s for workspaces with up to 500 data-links (SC-003)
**Constraints**: Read-only in this iteration; no cloud SDK on classpath (SC-006)
**Scale/Scope**: Single plugin change; no new plugin, no core module changes

---

## Constitution Check

| Principle | Status | Notes |
|-----------|--------|-------|
| I. Modular Architecture | ✅ PASS | Feature lives in `nf-tower`; no core module changes |
| II. Test-Driven Quality | ✅ PASS | Spock unit tests for client, handler, refactored path; reuses existing `Mock(TowerClient)` pattern |
| III. Dataflow Programming Model | ✅ PASS | NIO `InputStream`; no dataflow model changes |
| IV. Apache 2.0 License | ✅ PASS | All new files must include Apache 2.0 header |
| V. DCO Sign-off | ✅ PASS | All commits use `git commit -s` |
| VI. Semantic Versioning | ✅ PASS | VERSION bump and changelog entry both deferred to release time |
| VII. Groovy Idioms | ✅ PASS | `@CompileStatic`, follow existing `fs/` patterns |

No violations.

---

## Project Structure

### Documentation (this feature)

```text
specs/260422-seqera-datalinks-fs/
├── spec.md              ← feature spec (exists)
├── plan.md              ← this file
└── tasks.md             ← task checklist
```

### Source Code (nf-tower plugin)

```text
plugins/nf-tower/
└── src/     (VERSION and changelog.txt updated at release time, not in this feature)
    ├── main/io/seqera/tower/plugin/
    │   ├── fs/
    │   │   ├── ResourceTypeHandler.groovy         ← NEW (interface; list returns Iterable<Path>)
    │   │   ├── SeqeraFileSystemProvider.groovy    ← refactored (dispatch by handler; lazy filter iterator)
    │   │   ├── SeqeraFileSystem.groovy            ← refactored (handler registry; no dataset caches)
    │   │   ├── SeqeraPath.groovy                  ← refactored (trail segments, cachedAttributes, resolveWithAttributes)
    │   │   ├── SeqeraFileAttributes.groovy        ← refactored (isDir, size, lastMod, created, fileKey)
    │   │   ├── SeqeraPathFactory.groovy           ← unchanged
    │   │   ├── DatasetInputStream.groovy          ← unchanged
    │   │   └── handler/
    │   │         ├── DatasetsResourceHandler.groovy   ← NEW (extracted; owns dataset caches; parses @version)
    │   │         └── DataLinksResourceHandler.groovy  ← NEW
    │   ├── dataset/
    │   │   └── SeqeraDatasetClient.groovy         ← unchanged
    │   └── datalink/                              ← NEW package
    │         ├── SeqeraDataLinkClient.groovy      ← NEW (typed client; returns iterators and PagedDataLinkContent)
    │         └── PagedDataLinkContent.groovy      ← NEW (lazy pagination view over data-link content)
    └── test/io/seqera/tower/plugin/
        ├── fs/
        │   ├── SeqeraPathTest.groovy              ← extended (sub-path cases, cachedAttributes, trailing slash)
        │   ├── SeqeraFileSystemTest.groovy        ← extended (handler registry)
        │   ├── SeqeraFileSystemProviderTest.groovy ← extended (data-link dispatch specs)
        │   ├── ResourceTypeAbstractionTest.groovy ← NEW (architectural guard)
        │   └── handler/
        │         ├── DatasetsResourceHandlerTest.groovy   ← NEW (caches, attr short-circuit)
        │         └── DataLinksResourceHandlerTest.groovy  ← NEW (cache, credentialsId, paged listings)
        └── datalink/
              └── SeqeraDataLinkClientTest.groovy  ← NEW (pagination, endpoint URLs, error mapping)
```

**Structure decision**: Parallel `datalink/` package mirrors the existing `dataset/` package. Handlers live in `fs/handler/` so the generic NIO classes in `fs/` remain resource-type-agnostic. All wire DTOs are reused from `io.seqera.tower.model.*` — no plugin-local DTO classes. `PagedDataLinkContent` is a plugin-local service type (not a DTO) that wraps lazy pagination over `DataLinkItem` streams.

---

## Phase 0: Research Notes

### Tower-API DTOs (confirmed via `javap`)

All reused from `io.seqera:tower-api:1.121.0` (already on the classpath):

| DTO | Fields used here |
|---|---|
| `DataLinkDto` | `id: String`, `name: String`, `provider: DataLinkProvider`, `resourceRef: String`, `credentials: List<DataLinkCredentials>` |
| `DataLinkCredentials` | `id: String`, `name: String`, `provider: DataLinkProvider` |
| `DataLinkProvider` (enum) | `AWS`, `GOOGLE`, `AZURE`, `AZURE_ENTRA`, `AZURE_CLOUD`, `SEQERACOMPUTE`, `S3`. `toString()` returns the **lowercase** enum value (e.g. `"aws"`, `"google"`) — this is what appears in user-visible paths. |
| `DataLinksListResponse` | `dataLinks: List<DataLinkDto>`, `totalSize: Long` |
| `DataLinkContentResponse` | `originalPath: String`, `objects: List<DataLinkItem>`, `nextPageToken: String` |
| `DataLinkItem` | `type: DataLinkItemType`, `name: String`, `size: Long`, `mimeType: String` — no last-modified field |
| `DataLinkItemType` (enum) | `FOLDER`, `FILE` |
| `DataLinkDownloadUrlResponse` | `url: String` |

**Attribute consequence**: `DataLinkItem` does not expose a last-modified timestamp. `SeqeraFileAttributes.lastModifiedTime()` for data-link paths returns `FileTime.from(Instant.EPOCH)`. Spec assumption and FR-005 remain satisfied — we return a valid `FileTime`; the absence of real data is a Platform-API limitation.

### Platform endpoints (confirmed from OpenAPI)

| Operation | Endpoint | Notes |
|---|---|---|
| List data-links in workspace | `GET /data-links?workspaceId=<ws>&max=<n>&offset=<o>` | Offset pagination. `totalSize` = full count; `max=100` per page. Optional `&search=<name>` used by `getDataLink` for server-side pre-filter. |
| Browse root of a data-link | `GET /data-links/{id}/browse?workspaceId=<ws>` | Token pagination via `nextPageToken`. Optional `credentialsId`. |
| Browse a sub-path | `GET /data-links/{id}/browse/{path}?workspaceId=<ws>` | Same response and pagination as the root variant. The `{path}` segment preserves `/` as path separators. |
| Pre-signed download URL | `GET /data-links/{id}/generate-download-url?workspaceId=<ws>&filePath=<sub>` | Returns `DataLinkDownloadUrlResponse.url`. Optional `credentialsId`. |

`credentialsId` is taken from `DataLinkDto.credentials[0].id` when the list is non-empty, otherwise the query parameter is omitted.

### Signed-URL fetch

The signed URL is **not** a Seqera endpoint; it points at S3/GCS/Azure with auth baked into the query string. It must be fetched **without** the Seqera `Authorization` header (AWS SigV4 will reject unknown `Authorization` headers). Use a standalone `java.net.http.HttpClient` inside `DataLinksResourceHandler` for this fetch. Do **not** use `TowerClient.sendStreamingRequest`, which would add Seqera auth headers.

### SeqeraPath refactor shape

Replace the six typed fields (`org`, `workspace`, `resourceType`, `datasetName`, `version`, `relPath`) with:

- `org: String` (or null)
- `workspace: String` (or null)
- `resourceType: String` (or null)
- `trail: List<String>` (possibly empty) — the segments after `resourceType`
- `relPath: String` (for relative paths; mutually exclusive with absolute segments)
- `cachedAttributes: SeqeraFileAttributes` (nullable) — set only by handlers when this path is produced by a listing, so subsequent `readAttributes` calls skip the API

`trail` is opaque to `SeqeraPath` — handlers interpret it. Trail segments are stored verbatim (including any `@version` suffix for datasets); interpretation is the handler's responsibility. Concrete interpretations:

- **Dataset** (`resourceType = "datasets"`): `trail.size() == 0` → resource-type dir; `trail.size() == 1` → dataset file. The trail segment may carry an `@version` suffix (e.g. `samples@2`); `DatasetsResourceHandler.parseNameAndVersion` splits it internally.
- **Data-link** (`resourceType = "data-links"`): `trail.size() == 0` → resource-type dir; `trail.size() == 1` → provider dir; `trail.size() == 2` → data-link root dir; `trail.size() ≥ 3` → entry inside the data-link (directory or file, per `readAttributes`).

`depth()` becomes `3 + trail.size()` when `resourceType` is set, else the count of non-null identity fields.

`SeqeraPath` tolerates trailing slashes and accidental double-slashes in the URI (empty trail segments are filtered at parse time). `cachedAttributes` is ignored by `equals`/`hashCode`/`toString`/`toUri`/`resolve`/`getParent`; a new method `resolveWithAttributes(String segment, SeqeraFileAttributes attrs)` produces a child path carrying the given attrs.

### Existing tests to preserve

Running `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.*'` and `... --tests 'io.seqera.tower.plugin.dataset.*'` must continue to pass throughout the refactor. The dataset behavior does not change; only the class that implements it does.

---

## Phase 1: Design & Contracts

### `ResourceTypeHandler` interface

```groovy
interface ResourceTypeHandler {
    /** the depth-3 segment this handler owns, e.g. "datasets" or "data-links" */
    String getResourceType()

    /**
     * List entries at the given directory path. Caller has verified depth >= 3.
     * Returning Iterable lets implementations stream large listings without
     * materializing them in memory.
     */
    Iterable<Path> list(SeqeraPath dir) throws IOException

    /** return BasicFileAttributes for any path at depth >= 3 owned by this handler */
    SeqeraFileAttributes readAttributes(SeqeraPath path) throws IOException

    /** open a read stream for a leaf path; throw if the path is a directory */
    InputStream newInputStream(SeqeraPath path) throws IOException

    /** verify the path exists and modes are satisfiable; READ allowed, WRITE/EXECUTE rejected */
    void checkAccess(SeqeraPath path, AccessMode... modes) throws IOException
}
```

Handlers build each child path via `parent.resolveWithAttributes(segmentName, attrs)` so subsequent `readAttributes` calls short-circuit when the same path is used.

### `SeqeraDataLinkClient` contract

```groovy
class SeqeraDataLinkClient {
    SeqeraDataLinkClient(TowerClient towerClient)

    /**
     * Lazy iterator over every data-link in the workspace. Pages fetched on demand
     * via GET /data-links?workspaceId=<ws>&max=100&offset=<o>.
     */
    Iterator<DataLinkDto> listDataLinks(long workspaceId)

    /**
     * Server-side-filtered resolution of a single data-link by (provider, name).
     * Iterates /data-links with &search=<name>, short-circuits on first match;
     * result is @Memoized per (workspaceId, provider, name).
     * Throws NoSuchFileException if not found.
     */
    DataLinkDto getDataLink(long workspaceId, String provider, String name)

    /** Distinct provider identifiers present in the workspace (sorted). */
    Set<String> getDataLinkProviders(long workspaceId)

    /**
     * Lazy paginated view over /data-links/{id}/browse[/{path}].
     * The returned PagedDataLinkContent loads the first page eagerly and paginates
     * subsequent pages as its iterator advances.
     */
    PagedDataLinkContent getContent(String dataLinkId, String subPath, long workspaceId, String credentialsId = null)

    /** GET /data-links/{id}/generate-download-url?filePath=<sub>[&credentialsId=<c>] */
    DataLinkDownloadUrlResponse getDownloadUrl(String dataLinkId, String subPath, long workspaceId, String credentialsId = null)
}
```

All endpoints translate 401/403/404/5xx through the same `checkFsResponse` pattern used in `SeqeraDatasetClient`. The `credentialsId` parameter is forwarded as a query-string value when non-null; the handler sources it from `DataLinkDto.credentials[0].id`.

### `PagedDataLinkContent` contract

```groovy
class PagedDataLinkContent implements Iterable<DataLinkItem> {
    /** Page fetcher: fetch(null) -> first page; fetch(token) -> next page. */
    static interface PageFetcher {
        Map<String, Object> fetch(String nextPageToken) throws IOException
        // returns: {objects: List<DataLinkItem>, nextPageToken: String, originalPath: String (first page only)}
    }

    PagedDataLinkContent(String originalPath, List<DataLinkItem> firstPage, String firstPageNextToken, PageFetcher pageFetcher)

    String getOriginalPath()
    List<DataLinkItem> getFirstPage()        // eager, already loaded at construction
    boolean isEmpty()
    Iterator<DataLinkItem> iterator()        // yields first-page items, then paginates lazily
}
```

### `SeqeraFileSystem` handler registry

```groovy
class SeqeraFileSystem extends FileSystem {
    // existing org/workspace state (unchanged)
    private final Map<String, ResourceTypeHandler> handlers = new LinkedHashMap<>()

    void registerHandler(ResourceTypeHandler h) { handlers.put(h.resourceType, h) }
    ResourceTypeHandler getHandler(String type) { handlers.get(type) }
    Set<String> getResourceTypes() { Collections.unmodifiableSet(handlers.keySet()) }
}
```

`SeqeraFileSystemProvider.newFileSystem()` registers both handlers after constructing the filesystem:

```groovy
fs.registerHandler(new DatasetsResourceHandler(fs, new SeqeraDatasetClient(towerClient)))
fs.registerHandler(new DataLinksResourceHandler(fs, new SeqeraDataLinkClient(towerClient)))
```

### Dispatch in `SeqeraFileSystemProvider`

- Depth 0–2: handled directly (root/org/workspace) — uses `SeqeraFileSystem`'s org/workspace cache as before.
- Depth 2 listing: returns **`fs.getResourceTypes()`** as child paths (replaces the hard-coded `['datasets']`).
- Depth ≥ 3: dispatch to `fs.getHandler(sp.resourceType)`; unknown type → `NoSuchFileException("Unsupported resource type: ${sp.resourceType}")`.

### `SeqeraFileAttributes` refactor

Replace the `DatasetDto`-coupled constructor with two constructors:

```groovy
/** directory */
SeqeraFileAttributes(boolean isDir)
/** file with explicit metadata */
SeqeraFileAttributes(long size, Instant lastModified, Instant created, Object fileKey)
```

Internal fields become `(directory, size, lastModified, created, fileKey)`. `DatasetsResourceHandler` constructs the file variant from `DatasetDto`; `DataLinksResourceHandler` constructs it from `DataLinkItem`. The previous `SeqeraFileAttributes(DatasetDto)` and `SeqeraFileAttributes(boolean)` call sites are updated.

### `DataLinksResourceHandler` behaviors

| Path shape | Method | Implementation |
|---|---|---|
| `data-links/` (trail=[]) | `list` | `client.getDataLinkProviders(ws)` → distinct providers (sorted); emit child paths `data-links/<prov>` |
| `data-links/<prov>` (trail=[p]) | `list` | stream `client.listDataLinks(ws)`; collect names where provider matches; emit child paths; `NoSuchFileException` if none match |
| `data-links/<prov>/<name>` (trail=[p,n]) | `list` | `client.getDataLink(ws, p, n)` → `dl`; `client.getContent(dl.id, "", ws, credentialsIdOf(dl))` → wrap items as `Iterable<Path>` carrying cached `SeqeraFileAttributes` |
| `data-links/<prov>/<name>/<sub>/…` (trail ≥ 3) | `list` | same as above with `subPath = trail[2..].join('/')` |
| any depth ≥ 3 | `readAttributes` | short-circuit if `p.cachedAttributes` is set; else: data-link-root → directory; deeper → `getContent(id, sub, ws, credentialsIdOf(dl)).firstPage`; if one item matching the last segment with `type = FILE`, return file attrs (size from item); otherwise → directory |
| leaf file | `newInputStream` | `client.getDataLink(ws, p, n)` → `dl`; `client.getDownloadUrl(dl.id, sub, ws, credentialsIdOf(dl))`; open a plain JDK `HttpClient.send(..., BodyHandlers.ofInputStream())` against `response.url`; return body stream |

`credentialsIdOf(dl)` returns `dl.credentials[0].id` when non-empty, else `null` (query parameter omitted).

Provider segment canonicalization: the path segment is the `DataLinkProvider` enum's `toString()` — lowercase (e.g. `aws`, `google`, `azure`). A path with an unknown provider segment fails via `client.getDataLink(...)` → `NoSuchFileException`.

Listings populate cached attributes on each emitted `SeqeraPath` (via `parent.resolveWithAttributes(name, attrs)`) so a follow-up `readAttributes(child)` returns immediately with zero API calls. Attributes come directly from each `DataLinkItem`: file → `(size, Instant.EPOCH, Instant.EPOCH, item.name)`; folder → `SeqeraFileAttributes(true)`.

### Data-link identity resolution

`client.getDataLink(workspaceId, provider, name)` iterates `/data-links?search=<name>` (server-side pre-filter) and returns the first entry whose `provider.toString() == providerSegment`. Memoized via `@Memoized` per `(workspaceId, provider, name)` — repeated handler calls within a run hit the memoization cache. The handler does NOT maintain its own `Map<Long, List<DataLinkDto>>` cache — the client-level streaming iterator plus memoized-lookup replaces it.

---

## Phase 2: Deliverables

The detailed task list is in [tasks.md](tasks.md). Phases in execution order:

1. **Refactor (foundational)** — extract `ResourceTypeHandler`, `DatasetsResourceHandler`; generalize `SeqeraPath`, `SeqeraFileSystem`, `SeqeraFileSystemProvider`, `SeqeraFileAttributes`. Existing dataset tests must pass unchanged.
2. **Data-link API client** — implement `SeqeraDataLinkClient` with pagination and error mapping. Unit tests with `Mock(TowerClient)`.
3. **US1 — Read file inside a data-link** — implement `DataLinksResourceHandler.newInputStream` and register the handler.
4. **US2 — Browse hierarchy** — implement `list` and `readAttributes`; workspace listing enumerates handlers.
5. **US3 — Error paths** — explicit tests for unknown provider, unknown data-link, missing sub-path, 401/403 mapping.
6. **US4 — Extensibility validation** — architectural check that generic classes stay resource-type-agnostic; end-to-end spec exercising both handlers in one fs.
7. **Final verification** — full test suite green; no cloud SDK added to classpath. VERSION bump and changelog entry happen at release time, not in this feature.

Each task in `tasks.md` specifies exact file paths, exact code changes, exact test commands, and a commit step.
