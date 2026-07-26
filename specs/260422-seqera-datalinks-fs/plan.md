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
    │   │   ├── SeqeraFileSystemProvider.groovy    ← refactored (dispatch by handler; lazy filter iterator; one-shot DirectoryStream)
    │   │   ├── SeqeraFileSystem.groovy            ← refactored (holds TowerClient; cached getUserId; handler registry)
    │   │   ├── SeqeraPath.groovy                  ← refactored (trail segments, volatile cachedAttributes, resolveWithAttributes)
    │   │   ├── SeqeraFileAttributes.groovy        ← refactored (isDir, size, lastMod, created, fileKey)
    │   │   ├── SeqeraPathFactory.groovy           ← unchanged
    │   │   ├── DatasetInputStream.groovy          ← unchanged
    │   │   └── handler/
    │   │         ├── DatasetsResourceHandler.groovy   ← NEW (extracted; owns dataset caches; parses @version)
    │   │         └── DataLinksResourceHandler.groovy  ← NEW (parent-browse readAttributes; credentialsId forwarding; requireDataLink helper)
    │   ├── dataset/
    │   │   └── SeqeraDatasetClient.groovy         ← cleanup (only dataset endpoints; user/workspace lookup moved to SeqeraFileSystem)
    │   └── datalink/                              ← NEW package
    │         ├── SeqeraDataLinkClient.groovy      ← NEW (typed client; getDataLink uses combined keyword search; named static fetchers)
    │         └── PagedIterable.groovy             ← NEW (generic lazy pagination — eager first page, lazy subsequent)
    └── test/io/seqera/tower/plugin/
        ├── fs/
        │   ├── SeqeraPathTest.groovy              ← extended (sub-path cases, cachedAttributes, trailing slash)
        │   ├── SeqeraFileSystemTest.groovy        ← extended (TowerClient ctor, getUserId caching, handler registry)
        │   ├── SeqeraFileSystemProviderTest.groovy ← extended (data-link dispatch specs, DirectoryStream one-shot, cachedAttributes propagation)
        │   ├── ResourceTypeAbstractionTest.groovy ← NEW (architectural guard)
        │   └── handler/
        │         ├── DatasetsResourceHandlerTest.groovy   ← NEW (caches, attr short-circuit)
        │         └── DataLinksResourceHandlerTest.groovy  ← NEW (parent-browse, cache, credentialsId, paged listings)
        └── datalink/
              └── SeqeraDataLinkClientTest.groovy  ← NEW (pagination, getDataLink combined search, endpoint URLs, error mapping)
```

**Structure decision**: Parallel `datalink/` package mirrors the existing `dataset/` package. Handlers live in `fs/handler/` so the generic NIO classes in `fs/` remain resource-type-agnostic. User/workspace lookup lives on `SeqeraFileSystem` (shared infrastructure across resource types), not on a resource-type client. All wire DTOs are reused from `io.seqera.tower.model.*` — no plugin-local DTO classes. `PagedIterable<T>` is a plugin-local generic service type (not a DTO) that captures the eager-first-page + lazy-subsequent-page contract used by all paginated endpoints.

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
| User-id lookup | `GET /user-info` | Called from `SeqeraFileSystem.getUserId()`; result cached for the lifetime of the filesystem. |
| Workspaces for user | `GET /user/{userId}/workspaces` | Called from `SeqeraFileSystem.loadOrgWorkspaceCache()`. |
| List data-links in workspace | `GET /data-links?workspaceId=<ws>&max=<n>&offset=<o>` | Offset pagination. `totalSize` = full count; `max=100` per page. Optional `&search=<value>` keyword filter; `getDataLink` uses `&search=<name>+provider:<provider>` (URL-encoded) for a single server-side resolution. |
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
}
```

`checkAccess` is **not** on this interface — `SeqeraFileSystemProvider.checkAccess` rejects WRITE/EXECUTE upfront and delegates existence-check to `h.readAttributes(sp)`.

Handlers build each child path via `parent.resolveWithAttributes(segmentName, attrs)` so subsequent `readAttributes` calls short-circuit when the same path is used.

### `SeqeraDataLinkClient` contract

```groovy
class SeqeraDataLinkClient {
    SeqeraDataLinkClient(TowerClient towerClient)

    /**
     * Iterator over every data-link in the workspace, backed by PagedIterable.
     * The first page is fetched eagerly (early-fail at the call site); later pages
     * are fetched on demand. Endpoint: GET /data-links?workspaceId=<ws>&max=100&offset=<o>.
     */
    Iterator<DataLinkDto> listDataLinks(long workspaceId) throws IOException

    /**
     * Resolve a single data-link by (provider, name) using the Platform's combined
     * keyword search: `search=<name> provider:<provider>`. Server returns at most
     * one matching data-link; this method takes the first result and returns null
     * on miss. @Memoized per (ws, provider, name) — both successful resolutions
     * AND null misses are cached for the lifetime of the client.
     */
    DataLinkDto getDataLink(long workspaceId, String provider, String name) throws IOException

    /** Distinct provider identifiers present in the workspace (sorted, unmodifiable). @Memoized. */
    Set<String> getDataLinkProviders(long workspaceId)

    /**
     * Lazy paginated view over /data-links/{id}/browse[/{path}]. PagedIterable<DataLinkItem>
     * loads the first page eagerly and paginates subsequent pages on demand.
     * Optional &search=<value> for server-side prefix filter on entry names.
     */
    PagedIterable<DataLinkItem> getContent(String dataLinkId, String subPath, long workspaceId, String credentialsId = null, String search = null) throws IOException

    /** GET /data-links/{id}/generate-download-url?filePath=<sub>[&credentialsId=<c>] */
    DataLinkDownloadUrlResponse getDownloadUrl(String dataLinkId, String subPath, long workspaceId, String credentialsId = null)
}
```

All endpoints translate 401/403/404/5xx through the same `checkFsResponse` pattern used in `SeqeraDatasetClient`. The `credentialsId` parameter is forwarded as a query-string value when non-null; the handler sources it from `DataLinkDto.credentials[0].id`.

Pagination is delegated to two named static fetchers nested inside the client:
- `DataLinkListFetcher` — implements `PagedIterable.NextPageFetcher<DataLinkDto>`; offset-paginated; cursor state is `(offset, total)` instance fields.
- `DataLinkContentFetcher` — implements `PagedIterable.NextPageFetcher<DataLinkItem>`; token-paginated; cursor state is `nextToken` instance field.

### `PagedIterable<T>` contract

Generic lazy-pagination abstraction shared by all paginated endpoints. Eager first page (so `IOException` surfaces at the call site, not at the first `Iterator.hasNext()`); later pages on demand. Fetch failures during iteration wrap in `UncheckedIOException`.

```groovy
class PagedIterable<T> implements Iterable<T> {
    /** Stateful "give me the next page" callback. Cursor lives in the implementation. */
    static interface NextPageFetcher<T> {
        Page<T> fetch() throws IOException
    }

    static class Page<T> {
        final List<T> items
        final boolean isLast       // true → no more pages after this one
    }

    /** Eagerly fetch the first page; later pages on demand. */
    static <T> PagedIterable<T> start(NextPageFetcher<T> fetcher) throws IOException

    List<T> getFirstPage()         // eager, already loaded
    boolean isEmpty()
    Iterator<T> iterator()         // yields first-page items, then paginates lazily
}
```

### `SeqeraFileSystem` shape

Holds a `TowerClient` directly and exposes shared (across resource types) infrastructure:

```groovy
class SeqeraFileSystem extends FileSystem {
    SeqeraFileSystem(SeqeraFileSystemProvider provider, TowerClient towerClient)

    /** Cached for the lifetime of the FS — token doesn't change during a run. */
    long getUserId() throws IOException

    void loadOrgWorkspaceCache()
    Set<String> listOrgNames()
    List<String> listWorkspaceNames(String org)
    long resolveWorkspaceId(String org, String workspace) throws NoSuchFileException

    void registerHandler(ResourceTypeHandler h)
    ResourceTypeHandler getHandler(String type)
    Set<String> getResourceTypes()
}
```

`SeqeraFileSystemProvider.newFileSystem()` constructs the filesystem with the `TowerClient` and registers both handlers:

```groovy
fileSystem = new SeqeraFileSystem(this, towerClient)
fileSystem.registerHandler(new DatasetsResourceHandler(fileSystem, new SeqeraDatasetClient(towerClient)))
fileSystem.registerHandler(new DataLinksResourceHandler(fileSystem, new SeqeraDataLinkClient(towerClient)))
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
| `data-links/<prov>/<name>` (trail=[p,n]) | `list` | `requireDataLink(ws, p, n, dir)` → `dl`; `client.getContent(dl.id, "", ws, credentialsIdOf(dl))` → wrap items as `Iterable<Path>` carrying cached `SeqeraFileAttributes` |
| `data-links/<prov>/<name>/<sub>/…` (trail ≥ 3) | `list` | same as above with `subPath = trail[2..].join('/')` |
| `data-links/` or `data-links/<prov>` | `readAttributes` | trail=[] → directory; trail=[p] → validate via `client.getDataLinkProviders(ws)`; throw if unknown |
| `data-links/<prov>/<name>/<sub-path>` | `readAttributes` | short-circuit if `p.cachedAttributes` is set; else `requireDataLink(...)`; trail.size==2 → directory (data-link root); trail.size≥3 → **parent-browse**: list the parent directory and find the entry by name. Entry's `type` (FILE/FOLDER) is the authoritative signal; file → `(size, EPOCH, EPOCH, path)`; folder → directory marker; missing entry → `NoSuchFileException`. Server-side `&search=<lastSeg>` narrows the parent listing. |
| leaf file | `newInputStream` | `requireDataLink(...)`; `client.getDownloadUrl(dl.id, sub, ws, credentialsIdOf(dl))`; open a plain JDK `HttpClient.send(..., BodyHandlers.ofInputStream())` against `response.url`; return body stream |

`credentialsIdOf(dl)` returns `dl.credentials[0].id` when non-empty, else `null` (query parameter omitted).

`requireDataLink(ws, provider, name, pathForErrors)` calls `client.getDataLink(...)` and throws `NoSuchFileException` if the result is `null` — uniform error message across the three call sites.

Provider segment canonicalization: the path segment is the `DataLinkProvider` enum's `toString()` — lowercase (e.g. `aws`, `google`, `azure`). A path with an unknown provider segment fails via `client.getDataLink(...) == null` → `requireDataLink` → `NoSuchFileException`.

Listings populate cached attributes on each emitted `SeqeraPath` (via `parent.resolveWithAttributes(name, attrs)`) so a follow-up `readAttributes(child)` returns immediately with zero API calls. Attributes come directly from each `DataLinkItem`: file → `(size, Instant.EPOCH, Instant.EPOCH, item.name)`; folder → `SeqeraFileAttributes(true)`. The provider also writes back resolved attributes onto the path after a fresh `readAttributes`, so subsequent reads on the same instance also hit the cache.

### Data-link identity resolution

`client.getDataLink(workspaceId, provider, name)` issues a single combined keyword search (`&search=<name>+provider:<provider>`, URL-encoded) — the Platform returns at most the matching data-link, and the method returns the first result or `null`. `@Memoized` per `(workspaceId, provider, name)` — both successful resolutions AND `null` misses are cached for the lifetime of the client. The handler does NOT maintain its own data-link cache; the client-level memoization replaces it. Why no client-side iterate-and-filter: the server-side `provider:<x>` keyword filter eliminates the need.

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
