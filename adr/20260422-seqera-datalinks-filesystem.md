# NIO Filesystem Support for Seqera Platform Data-Links

- Authors: Jorge Ejarque
- Status: draft
- Date: 2026-04-22
- Tags: nio, filesystem, seqera, data-links, nf-tower

Technical Story: Extend the `seqera://` NIO filesystem (introduced by [20260310-seqera-dataset-filesystem](20260310-seqera-dataset-filesystem.md)) to address files and directories inside Seqera Platform data-links without requiring cloud-provider credentials or SDK integration.

## Summary

Add a second resource type (`data-links`) to the existing `seqera://` filesystem in the `nf-tower` plugin. Paths of the form `seqera://<org>/<ws>/data-links/<provider>/<name>/<sub-path>` resolve to files and directories inside a Platform-managed data-link. Listings and attribute queries are served by the Platform's `/data-links/{id}/content` endpoint; byte reads go through pre-signed URLs returned by `/data-links/{id}/download`. Only the Seqera access token is required — no AWS/GCP/Azure credentials, no cloud SDK dependency.

As part of this change, the existing dataset-specific logic in `SeqeraFileSystemProvider`, `SeqeraFileSystem`, and `SeqeraPath` is extracted into a `ResourceTypeHandler` abstraction, so the two resource types coexist behind a common contract.

## Problem Statement

The dataset filesystem ships `seqera://` URI support for Platform datasets, but datasets are only one of several file-like resources users manage on the Platform. Data-links are the most common — they reference a cloud bucket or prefix (S3, GCS, Azure Blob) with potentially large, nested content. Today, a pipeline that needs to read a file inside a data-link must:

1. Look up the data-link's underlying URI outside Nextflow.
2. Configure cloud credentials in the compute environment (AWS access keys, GCP service account, Azure SAS, etc.).
3. Reference the object by its cloud URI.

This is friction the Platform already solves: data-links are scoped, ACL-controlled entities, and the Platform knows how to broker access to their content. A `seqera://` URI for a path inside a data-link would let pipelines consume Platform-managed data with only the Seqera access token — no cloud SDK, no credential sprawl.

## Goals or Decision Drivers

- Native `seqera://` access to files and directories inside Platform data-links, at arbitrary depth.
- Zero cloud-provider credential configuration — the Seqera access token is the only auth surface.
- No new runtime dependency on cloud SDKs (`aws-sdk`, `google-cloud-storage`, `azure-*`).
- Reuse of existing nf-tower plugin infrastructure — `TowerClient` for HTTP + auth + retry, tower-api DTOs for wire types.
- Introduce a `ResourceTypeHandler` abstraction so the dataset and data-link behaviors share one filesystem without leaking into each other.
- Preserve Platform-side access control for listings and metadata (not just reads).

## Non-goals

- Write operations to data-links (upload). The Platform's `POST /data-links/{id}/multipart-upload` is a future hook; it is not implemented in this iteration.
- Data-link management (create/update/delete the data-link entity itself).
- Transparent pre-signed URL renewal when a URL expires mid-stream — failures surface as `IOException` and Nextflow task retry handles them.
- Browse-result caching within a run.
- Fusion integration — Fusion has its own data-link access path.

## Considered Options

### Option 1: Platform-brokered credentials + cloud SDK delegation

For each read, call the Platform to obtain short-lived AWS/GCP/Azure credentials scoped to the data-link, then use the existing `nf-amazon` / `nf-google` / `nf-azure` providers for the actual I/O.

- Good, because cloud providers handle streaming, range reads, multi-part efficiently.
- Bad, because it requires `nf-tower` to depend on (or coordinate with) three cloud plugins.
- Bad, because credential plumbing across plugin boundaries is complex — each cloud plugin has its own credential object model.
- Bad, because it adds failure modes around credential refresh windows crossing long reads.

### Option 2: Pre-signed URL + direct HTTPS fetch

Call the Platform's `GET /data-links/{id}/download?path=<sub>` endpoint to obtain a pre-signed URL; stream bytes through the existing `TowerClient` / `HxClient` HTTPS path.

- Good, because there is no cloud SDK dependency — all I/O is generic HTTPS.
- Good, because the Platform is the only credential surface (user token goes in, signed URL comes out; credentials never cross our process boundary as a distinct object).
- Good, because it uniformly supports every provider the Platform supports — now and in the future — with no per-provider code.
- Good, because `TowerClient.sendStreamingRequest()` already exists and has the retry/backoff semantics we want.
- Bad, because pre-signed URLs have time windows; a very long read can outlive its URL. Acceptable: Nextflow task retry handles the failure.
- Bad, because range reads / multi-part reads are not implemented in this iteration. Acceptable: datasets are already single-shot reads and the pattern matches.

### Option 3: Proxy all bytes through the Platform

Route all reads through a Platform endpoint that streams content back to the client from the underlying cloud.

- Good, because the Platform sees and can log every byte.
- Bad, because it imposes Platform bandwidth/egress cost on every pipeline byte.
- Bad, because no such primary endpoint is offered — `/download` returns a URL, not bytes.

## Pros and Cons of the Options

See above.

## Solution or decision outcome

Option 2 — pre-signed URL + direct HTTPS fetch. All data-link byte I/O goes through the Platform's `/download` endpoint and `TowerClient.sendStreamingRequest()`. The plugin never touches a cloud SDK and never holds a long-lived cloud credential.

Extend the `fs/` package with a real `ResourceTypeHandler` abstraction. Extract the existing dataset logic into a `DatasetsResourceHandler`. Add `DataLinksResourceHandler` as the second implementation.

## Rationale & discussion

### Path Hierarchy

The `seqera://` path gains a second resource-type branch:

```
seqera://                                                     → ROOT (directory, depth 0)
  └── <org>/                                                  → ORGANIZATION (directory, depth 1)
        └── <workspace>/                                      → WORKSPACE (directory, depth 2)
              ├── datasets/                                   → RESOURCE TYPE (directory, depth 3)
              │     └── <name>[@<version>]                    → DATASET (file, depth 4)
              └── data-links/                                 → RESOURCE TYPE (directory, depth 3)
                    └── <provider>/                           → PROVIDER (directory, depth 4)
                          └── <name>/                         → DATA-LINK (directory, depth 5)
                                └── <sub>/…/<file-or-dir>     → CONTENT (directory or file, depth 6+)
```

Three structural differences from datasets:

1. **Two identity segments** (`<provider>/<name>`) instead of one (`<name>`). Provider disambiguation is required because a workspace can host two data-links with the same name on different clouds.
2. **Arbitrary sub-path depth** below the data-link root. Each segment is a folder or file inside the underlying bucket.
3. **No version pinning** — data-link content is not versioned by the Platform. Content is always "current".

`ResourceTypeHandler.getIdentitySegmentCount()` encodes the difference: 1 for datasets, 2 for data-links. `SeqeraPath` treats everything after the identity segments as the handler-owned sub-path.

### Component Structure

```
plugins/nf-tower/src/main/io/seqera/tower/plugin/
├── fs/                                   ← generic NIO layer (refactored)
│   ├── SeqeraFileSystemProvider          ← dispatches by resourceType to handler
│   ├── SeqeraFileSystem                  ← org/ws cache + handler registry
│   ├── SeqeraPath                        ← generic segment list (identity + sub-path)
│   ├── SeqeraFileAttributes              ← plain (isDir, size, lastModified) holder
│   ├── SeqeraPathFactory                 ← unchanged
│   ├── DatasetInputStream                ← unchanged
│   ├── ResourceTypeHandler               ← NEW interface
│   └── handler/
│         ├── DatasetsResourceHandler     ← NEW — dataset logic extracted here
│         └── DataLinksResourceHandler    ← NEW
├── dataset/
│   └── SeqeraDatasetClient               ← unchanged
└── datalink/                             ← NEW
      └── SeqeraDataLinkClient            ← typed client over TowerClient
                                             returns io.seqera.tower.model.* directly
```

No plugin-local DTO classes are introduced. `DataLinkDto`, `DataLinkContentResponse`, `DataLinkItem`, `DataLinkDownloadUrlResponse`, `DataLinkProvider` and related types are reused from `io.seqera:tower-api:1.121.0`.

### `ResourceTypeHandler` contract

```
interface ResourceTypeHandler {
    String getResourceType()                             // "datasets" / "data-links"
    int getIdentitySegmentCount()                        // 1 / 2
    List<Path> list(SeqeraPath dir) throws IOException
    SeqeraFileAttributes readAttributes(SeqeraPath p) throws IOException
    InputStream newInputStream(SeqeraPath p) throws IOException
    void checkAccess(SeqeraPath p, AccessMode... modes) throws IOException
}
```

`SeqeraFileSystemProvider` owns dispatch at depth ≥ 3. Depth 0–2 (root/org/workspace) remains in `SeqeraFileSystem`, shared across all handlers. At depth 3 (the workspace listing returns the resource-type children), the handler registry is enumerated — `datasets` and `data-links` are the two entries today, added automatically by the provider at `newFileSystem()` time.

### API Usage Summary (Data-Links)

| NIO operation                                                  | Platform endpoint                          | Notes                                                                           |
| -------------------------------------------------------------- | ------------------------------------------ | ------------------------------------------------------------------------------- |
| list data-links in workspace (resolution + depth-4/5 listings) | `GET /data-links?workspaceId=X`            | cached per workspace; pagination exhausted                                      |
| `newDirectoryStream(dir)` at depth ≥ 5                         | `GET /data-links/{id}/content?path=<sub>`  | items array → child paths                                                       |
| `readAttributes(path)` inside a data-link                      | `GET /data-links/{id}/content?path=<sub>`  | single targeted call; file vs folder from response                              |
| `newInputStream(file)`                                         | `GET /data-links/{id}/download?path=<sub>` | parse `DataLinkDownloadUrlResponse.url`; fetch with plain JDK `HttpClient` (no Seqera auth header — the URL is signed for the cloud backend) |

### Key Design Decisions

1. **TowerClient delegation for Platform calls**: `SeqeraDataLinkClient` routes all Seqera API calls (list, content, download-URL) through `TowerClient.sendApiRequest()`, sharing authentication state with the dataset client. The pre-signed URL itself is fetched directly with a plain JDK `HttpClient` — no Seqera headers are sent to the cloud backend.

2. **Pre-signed URLs, not credential brokering**: the Platform returns a URL that already has the auth embedded. No AWS/GCP/Azure SDK is imported; no credential object crosses the plugin boundary. This is the single biggest simplification relative to a "get creds, hand to cloud plugin" approach.

3. **No per-stream URL renewal**: if a signed URL expires mid-read, the HTTP connection errors and the `InputStream` surfaces an `IOException`. Nextflow task retry handles the failure as it does for any other transient read failure. The plugin does not implement transparent re-issuance.

4. **Provider disambiguation in the path**: the data-link identity is `(workspace, provider, name)` on the Platform side. The path segment layout mirrors this to avoid ambiguity when names collide across providers.

5. **Reuse tower-api DTOs**: every wire type is an `io.seqera.tower.model.*` class already on the plugin's classpath via `tower-api:1.121.0`. No parallel plugin-local DTOs.

6. **Handler registry at construction, not via PF4J**: handlers are instantiated in `SeqeraFileSystemProvider.newFileSystem()`. Adding a third resource type is a code change to this plugin, identical in shape to the dataset/data-link pair. No extension-point protocol is introduced — YAGNI.

7. **`readAttributes` is single-target**: because `GET /data-links/{id}/content?path=<sub>` accepts both directory and file paths, a file-level `readAttributes` is one API call — not a parent browse plus filter. No N+1 problem; no browse cache needed.

8. **Read-only stance preserved**: `SeqeraFileSystem.isReadOnly()` remains `true`. Write operations on data-links raise `UnsupportedOperationException`. The `/data-links/{id}/multipart-upload` endpoint is a future extension point.

9. **Minimal cache**: only the workspace-level data-link list is cached. No browse-result cache. No URL cache. Simpler to reason about; consistent with the dataset handler's cache shape.

### Refactor Delivered by This Change

Adding a second resource type requires a shared abstraction in the `fs/` package so the two behaviors do not collide:

- The `ResourceTypeHandler` interface is introduced.
- All dataset-specific logic previously inlined in `SeqeraFileSystemProvider`, `SeqeraFileSystem`, and `SeqeraPath` is moved to a new `DatasetsResourceHandler`.
- `DataLinksResourceHandler` is added alongside it, implementing the same interface.
- The generic classes (`SeqeraFileSystemProvider`, `SeqeraFileSystem`, `SeqeraPath`) become resource-type-agnostic for depth ≥ 3 — they dispatch to handlers and carry no knowledge of either resource's semantics.

The existing dataset test suite continues to pass unchanged; every dataset code path is routed through `DatasetsResourceHandler` without behavioral change.

### Limitations

- **No write support for data-links in this iteration.** Upload paths must continue to use Fusion or direct cloud-SDK access until a follow-up adds the multipart-upload handler.
- **Signed URL expiration is not handled transparently.** Very long reads may outlive the URL's validity window.
- **No transparent pagination of data-link entries inside a single directory.** The browse API's pagination (if any) must be exhausted; for very large directories this may be slower than direct cloud listings.
- **Single endpoint per JVM** (unchanged from dataset ADR): concurrent access to multiple Platform endpoints in one JVM is not supported.

## Links

- [Spec](../specs/260422-seqera-datalinks-fs/spec.md)
- Extends [20260310-seqera-dataset-filesystem](20260310-seqera-dataset-filesystem.md)

## More information

- [Seqera Platform OpenAPI spec](https://cloud.seqera.io/openapi/seqera-api-latest.yml) — `/data-links` endpoints.
- [What is an ADR and why should you use them](https://github.com/thomvaill/log4brains/tree/master#-what-is-an-adr-and-why-should-you-use-them)
