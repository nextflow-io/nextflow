# Implementation Plan: Seqera NIO Filesystem Support for Platform Data-Links

**Branch**: `260422-seqera-datalinks-fs` | **Date**: 2026-04-23 | **Spec**: [spec.md](spec.md) | **ADR**: [20260422-seqera-datalinks-filesystem](../../adr/20260422-seqera-datalinks-filesystem.md)

---

## Summary

Extend the `nf-tower` plugin's `seqera://` NIO filesystem (shipped in [260310-seqera-dataset-fs](../260310-seqera-dataset-fs/spec.md)) with a second resource type, `data-links`. Paths of the form `seqera://<org>/<ws>/data-links/<provider>/<name>/<sub-path>` resolve to entries inside Platform-managed data-links (S3/GCS/Azure buckets or prefixes). Listings and attribute queries hit the Platform's `/data-links/{id}/content` endpoint; byte reads go through a pre-signed URL returned by `/data-links/{id}/download` and fetched with a plain JDK HTTP client — no cloud SDK dependency.

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
    │   │   ├── ResourceTypeHandler.groovy         ← NEW (interface)
    │   │   ├── SeqeraFileSystemProvider.groovy    ← refactored (dispatch by handler)
    │   │   ├── SeqeraFileSystem.groovy            ← refactored (handler registry)
    │   │   ├── SeqeraPath.groovy                  ← refactored (generic trail segments)
    │   │   ├── SeqeraFileAttributes.groovy        ← refactored (isDir, size, lastMod)
    │   │   ├── SeqeraPathFactory.groovy           ← unchanged
    │   │   ├── DatasetInputStream.groovy          ← unchanged
    │   │   └── handler/
    │   │         ├── DatasetsResourceHandler.groovy   ← NEW (extracted)
    │   │         └── DataLinksResourceHandler.groovy  ← NEW
    │   ├── dataset/
    │   │   └── SeqeraDatasetClient.groovy         ← unchanged
    │   └── datalink/                              ← NEW package
    │         └── SeqeraDataLinkClient.groovy      ← NEW
    └── test/io/seqera/tower/plugin/
        ├── fs/
        │   ├── SeqeraPathTest.groovy              ← extended (sub-path cases)
        │   ├── SeqeraFileSystemTest.groovy        ← extended (handler registry)
        │   ├── SeqeraFileSystemProviderTest.groovy ← extended (data-link dispatch specs)
        │   └── handler/
        │         ├── DatasetsResourceHandlerTest.groovy   ← NEW
        │         └── DataLinksResourceHandlerTest.groovy  ← NEW
        └── datalink/
              └── SeqeraDataLinkClientTest.groovy  ← NEW
```

**Structure decision**: Parallel `datalink/` package mirrors the existing `dataset/` package. Handlers live in `fs/handler/` so the generic NIO classes in `fs/` remain resource-type-agnostic. All wire DTOs are reused from `io.seqera.tower.model.*` — no plugin-local DTO classes.

---

## Phase 0: Research Notes

### Tower-API DTOs (confirmed via `javap`)

All reused from `io.seqera:tower-api:1.121.0` (already on the classpath):

| DTO | Fields used here |
|---|---|
| `DataLinkDto` | `id: String`, `name: String`, `provider: DataLinkProvider`, `resourceRef: String` |
| `DataLinkProvider` (enum) | `AWS`, `GOOGLE`, `AZURE`, `AZURE_ENTRA`, `AZURE_CLOUD`, `SEQERACOMPUTE`, `S3` — exposes a `String value` via `toString()` |
| `DataLinksListResponse` | `dataLinks: List<DataLinkDto>`, `totalSize: Long` |
| `DataLinkContentResponse` | `originalPath: String`, `objects: List<DataLinkItem>`, `nextPageToken: String` |
| `DataLinkItem` | `type: DataLinkItemType`, `name: String`, `size: Long`, `mimeType: String` — no last-modified field |
| `DataLinkItemType` (enum) | `FOLDER`, `FILE` |
| `DataLinkDownloadUrlResponse` | `url: String` |

**Attribute consequence**: `DataLinkItem` does not expose a last-modified timestamp. `SeqeraFileAttributes.lastModifiedTime()` for data-link paths returns `FileTime.from(Instant.EPOCH)`. Spec assumption 2 and FR-005 remain satisfied — we return a valid `FileTime`; the absence of real data is a Platform-API limitation.

### Platform endpoints (confirmed from OpenAPI)

- `GET /data-links?workspaceId=<wsId>&max=<n>&offset=<n>` → `DataLinksListResponse`. Pagination: `totalSize` is the full count; keep fetching with offset until sum of received `dataLinks` equals `totalSize`. Default `max` is the server default; plugin uses `max=100` per page.
- `GET /data-links/{id}/content?workspaceId=<wsId>&path=<sub>&nextPageToken=<tok>` → `DataLinkContentResponse`. Works for directory and file paths. Pagination: follow `nextPageToken` until null/empty.
- `GET /data-links/{id}/download?workspaceId=<wsId>&path=<sub>` → `DataLinkDownloadUrlResponse` with a cloud-signed URL.

### Signed-URL fetch

The signed URL is **not** a Seqera endpoint; it points at S3/GCS/Azure with auth baked into the query string. It must be fetched **without** the Seqera `Authorization` header (AWS SigV4 will reject unknown `Authorization` headers). Use a standalone `java.net.http.HttpClient` inside `DataLinksResourceHandler` for this fetch. Do **not** use `TowerClient.sendStreamingRequest`, which adds Seqera auth headers.

### SeqeraPath refactor shape

Replace the six typed fields (`org`, `workspace`, `resourceType`, `datasetName`, `version`, `relPath`) with:

- `org: String` (or null)
- `workspace: String` (or null)
- `resourceType: String` (or null)
- `trail: List<String>` (possibly empty) — the segments after `resourceType`
- `relPath: String` (for relative paths; mutually exclusive with absolute segments)

The `trail` is opaque to `SeqeraPath` — handlers interpret it. Concrete interpretations:

- **Dataset** (`resourceType = "datasets"`): `trail.size() == 0` → resource-type dir; `trail.size() == 1` → dataset file, with optional `@version` suffix on the single element; `trail.size() > 1` → invalid.
- **Data-link** (`resourceType = "data-links"`): `trail.size() == 0` → resource-type dir; `trail.size() == 1` → provider dir; `trail.size() == 2` → data-link root dir; `trail.size() ≥ 3` → entry inside the data-link (directory or file, per `readAttributes`).

`depth()` becomes `3 + trail.size()` when `resourceType` is set, else the count of non-null identity fields.

### Existing tests to preserve

Running `./gradlew :plugins:nf-tower:test --tests 'io.seqera.tower.plugin.fs.*'` and `... --tests 'io.seqera.tower.plugin.dataset.*'` must continue to pass throughout the refactor. The dataset behavior does not change; only the class that implements it does.

---

## Phase 1: Design & Contracts

### `ResourceTypeHandler` interface

```groovy
interface ResourceTypeHandler {
    /** the depth-3 segment this handler owns, e.g. "datasets" or "data-links" */
    String getResourceType()

    /** list entries at the given directory path owned by this handler; caller verified depth >= 3 and isDirectory */
    List<Path> list(SeqeraPath dir) throws IOException

    /** return BasicFileAttributes for any path at depth >= 3 owned by this handler */
    SeqeraFileAttributes readAttributes(SeqeraPath path) throws IOException

    /** open a read stream for a leaf path; throw if the path is a directory */
    InputStream newInputStream(SeqeraPath path) throws IOException

    /** verify the path exists and modes are satisfiable; READ allowed, WRITE/EXECUTE rejected */
    void checkAccess(SeqeraPath path, AccessMode... modes) throws IOException
}
```

### `SeqeraDataLinkClient` contract

```groovy
class SeqeraDataLinkClient {
    SeqeraDataLinkClient(TowerClient towerClient)

    /** exhaust pagination; return all data-links in the workspace */
    List<DataLinkDto> listDataLinks(long workspaceId)

    /** GET /data-links/{id}/content?path=<sub> — exhausts nextPageToken pagination */
    DataLinkContentResponse getContent(String dataLinkId, String subPath, long workspaceId)

    /** GET /data-links/{id}/download?path=<sub> */
    DataLinkDownloadUrlResponse getDownloadUrl(String dataLinkId, String subPath, long workspaceId)
}
```

All three translate 401/403/404/5xx through the same `checkFsResponse` pattern used in `SeqeraDatasetClient`.

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
| `data-links/` (trail=[]) | `list` | enumerate distinct `DataLinkDto.provider` (via `toString()` / enum value) in cached list; return paths `data-links/<prov>` |
| `data-links/<prov>` (trail=[p]) | `list` | filter cached list where provider matches; return paths `data-links/<prov>/<name>` |
| `data-links/<prov>/<name>` (trail=[p,n]) | `list` | resolve to `dataLinkId`; call `getContent(id, "", wsId)`; map `objects` to child paths |
| `data-links/<prov>/<name>/<sub>/…` (trail ≥ 3) | `list` | call `getContent(id, "<sub>/…", wsId)`; map `objects` |
| any depth ≥ 3 | `readAttributes` | data-link-root path → directory; below that → `getContent(id, sub, ws)`; if response `objects` has one item matching the last segment with `type = FILE`, return file attrs (size from item); otherwise → directory |
| leaf file | `newInputStream` | `getDownloadUrl(id, sub, ws)`; open plain `HttpClient.send(..., BodyHandlers.ofInputStream())` against `response.url`; return body stream |

Provider segment canonicalization: the path segment is the `DataLinkProvider` enum's `toString()`. A path with an unknown provider segment maps to `NoSuchFileException`.

### Data-link identity resolution

`DataLinksResourceHandler.resolveDataLinkId(provider, name, workspaceId)`:

1. Ensure the workspace's data-link list is loaded (cached `Map<Long, List<DataLinkDto>>` inside the handler).
2. Find `DataLinkDto` where `provider.toString() == providerSegment && name == nameSegment`.
3. Return `id`; throw `NoSuchFileException` with a clear message if not found.

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
