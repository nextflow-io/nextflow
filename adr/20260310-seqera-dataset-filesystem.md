# NIO Filesystem for Seqera Platform Datasets

- Authors: Jorge Ejarque
- Status: draft
- Date: 2026-03-10
- Tags: nio, filesystem, seqera, datasets, nf-tower

Technical Story: Enable Nextflow pipelines to read Seqera Platform datasets as ordinary file paths using `seqera://` URIs.

## Summary

Add a Java NIO `FileSystemProvider` to the `nf-tower` plugin that registers the `seqera://` scheme, allowing pipelines to reference Seqera Platform datasets (CSV/TSV) as standard file paths without manual download steps. The implementation reuses the existing `TowerClient` for all HTTP communication, inheriting authentication and retry behaviour.

## Problem Statement

Nextflow users managing datasets on the Seqera Platform must currently download dataset files manually or through custom scripts before referencing them in pipelines. There is no native integration between Nextflow's file abstraction and the Seqera Platform dataset API. This creates friction in workflows where datasets are the primary input and forces users to handle authentication, versioning, and file staging outside the pipeline definition.

## Goals or Decision Drivers

- Transparent access to Seqera Platform datasets using standard Nextflow file path syntax
- Reuse of existing nf-tower plugin infrastructure (authentication, HTTP client, retry/backoff)
- Hierarchical path browsing matching the platform's org/workspace/dataset structure
- Extensible architecture that can support future Seqera-managed resource types (e.g. data-links)
- No new plugin or module — feature lives within nf-tower

## Non-goals

- Streaming large datasets — the Platform API does not support streaming; content is fully buffered on download
- Implementing resource types beyond `datasets` — only the extensible architecture is required
- Local caching across pipeline runs — Nextflow's standard task staging handles caching
- Dataset management operations (delete, rename) — the filesystem is read-only in the initial implementation

## Considered Options

### Option 1: Standalone plugin with dedicated HTTP client

A new `nf-seqera-fs` plugin with its own HTTP client configuration and authentication setup.

- Good, because it isolates the filesystem code from the nf-tower plugin
- Bad, because it duplicates authentication configuration and HTTP client setup
- Bad, because two separate HTTP clients sharing a refresh token would corrupt each other's auth state

### Option 2: NIO filesystem within nf-tower using TowerClient delegation

Add the filesystem to nf-tower, delegating all HTTP through the existing `TowerClient` singleton via a typed `SeqeraDatasetClient` wrapper.

- Good, because it shares authentication and token refresh with TowerClient
- Good, because it reuses existing retry/backoff configuration
- Good, because no new dependencies are needed

### Option 3: Direct HxClient usage within nf-tower

Add the filesystem to nf-tower but use `HxClient` directly rather than going through TowerClient.

- Good, because it gives full control over request construction
- Bad, because exposing HxClient internals couples the filesystem to implementation details
- Bad, because token refresh coordination with TowerClient becomes manual

## Solution or decision outcome

Option 2 — NIO filesystem within nf-tower using TowerClient delegation. All HTTP calls go through `TowerClient.sendApiRequest()`, ensuring a single point of authentication and retry logic.

## Rationale & discussion

### Path Hierarchy

The `seqera://` path encodes the Platform's organizational structure directly:

```
seqera://                                        → ROOT (directory, depth 0)
  └── <org>/                                     → ORGANIZATION (directory, depth 1)
        └── <workspace>/                         → WORKSPACE (directory, depth 2)
              └── datasets/                      → RESOURCE TYPE (directory, depth 3)
                    └── <name>[@<version>]        → DATASET (file, depth 4)
```

Each level is a directory except the leaf dataset, which is a file. Version pinning uses an `@version` suffix on the dataset name segment (e.g. `seqera://acme/research/datasets/samples@2`). Without it, the latest non-disabled version is resolved.

### Name-to-ID Resolution

The path uses human-readable names but the Platform API requires numeric IDs. Resolution is built from two API calls at filesystem initialization:

1. `GET /user-info` → obtain `userId`
2. `GET /user/{userId}/workspaces` → returns all accessible org/workspace pairs

This single source provides both directory listing content and name→ID mapping. Results are cached in `SeqeraFileSystem` with invalidation on write operations. `GET /orgs` is intentionally not used as it returns all platform orgs, not scoped to user membership.

### Component Structure

```
plugins/nf-tower/src/main/io/seqera/tower/plugin/
├── fs/                             ← NIO layer
│   ├── SeqeraFileSystemProvider    ← FileSystemProvider (scheme: "seqera")
│   ├── SeqeraFileSystem            ← FileSystem with org/workspace/dataset caches
│   ├── SeqeraPath                  ← Path implementation (depth 0–4)
│   ├── SeqeraFileAttributes        ← BasicFileAttributes
│   ├── SeqeraPathFactory           ← PF4J FileSystemPathFactory extension
│   └── DatasetInputStream          ← SeekableByteChannel over InputStream
├── dataset/                        ← API client layer
│   ├── SeqeraDatasetClient         ← Typed HTTP client wrapping TowerClient
│   ├── DatasetDto                  ← Dataset API response model
│   ├── DatasetVersionDto           ← Version API response model
│   ├── OrgAndWorkspaceDto          ← Org/workspace list model
│   └── WorkspaceOrgDto             ← Workspace/org mapping model
└── resources/META-INF/services/
    └── java.nio.file.spi.FileSystemProvider
```

### Key Design Decisions

1. **TowerClient delegation**: `SeqeraDatasetClient` delegates all HTTP through `TowerFactory.client()` → `TowerClient.sendApiRequest()`. This ensures shared authentication state and avoids the token refresh corruption that would occur with separate HTTP client instances.

2. **One filesystem per JVM**: `SeqeraFileSystemProvider` maintains a single `SeqeraFileSystem` keyed by scheme. This matches the `TowerClient` singleton-per-session pattern.

3. **Read-only initial scope**: The filesystem reports `isReadOnly()=true`. Write support (dataset upload via multipart POST) is deferred to a future iteration.

4. **Download filename constraint**: The Platform API's download endpoint (`GET /datasets/{id}/v/{version}/n/{fileName}`) requires the exact filename from upload time. The implementation always resolves `DatasetVersionDto.fileName` from `GET /datasets/{id}/versions` before constructing the download URL.

5. **Extensible resource types**: The path hierarchy reserves depth 3 for a resource type segment (currently only `datasets`). Adding support for data-links or other resource types requires only a new handler at the directory listing and I/O layers, with no changes to path resolution or authentication.

6. **Thread safety**: `SeqeraFileSystem` cache methods and `SeqeraFileSystemProvider` lifecycle methods are `synchronized`. The filesystem map uses `LinkedHashMap` with external synchronization rather than `ConcurrentHashMap`, matching the low-contention access pattern.

### Limitations

- **No size metadata**: `SeqeraFileAttributes.size()` returns 0 for all paths because the Platform API does not expose content length in dataset metadata.
- **Single endpoint per JVM**: The filesystem key is scheme-only; concurrent access to different Platform endpoints in the same JVM is not supported.

### Streaming Downloads

Dataset downloads use `TowerClient.sendStreamingRequest()` which calls `HxClient.sendAsStream()` — the response body is returned as an `InputStream` streamed directly from the HTTP connection. This avoids the triple-buffering problem (`String` → `getBytes()` → `ByteArrayInputStream`) that would otherwise consume ~40 MB heap per 10 MB dataset. The `HxClient.sendAsStream()` method goes through the same `sendWithRetry()` path as `sendAsString()`, so retry logic and token refresh are preserved.

## Links

- [Spec](../specs/260310-seqera-dataset-fs/spec.md)
- [Implementation plan](../specs/260310-seqera-dataset-fs/plan.md)
- [Data model](../specs/260310-seqera-dataset-fs/data-model.md)
