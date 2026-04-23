# Feature Specification: Seqera NIO Filesystem Support for Platform Data-Links

**Feature Branch**: `260422-seqera-datalinks-fs`
**Created**: 2026-04-22
**Status**: Draft
**Depends on**: [260310-seqera-dataset-fs](../260310-seqera-dataset-fs/spec.md)
**Input**: User description: "I want to extend the seqera NIO filesystem to include the seqera platform data-links using this url template `seqera://<org>/<workspace>/data-links/...`. The seqera platform API is https://cloud.seqera.io/openapi/seqera-api-latest.yml"

## Clarifications

### Session 2026-04-22

- Q: Scope of data-link support — list-only, full Platform-driven traversal, or hybrid (Platform-driven listing + cloud-driven I/O)? → A: Hybrid. Listing and attributes go through the Seqera Platform; byte-level I/O is delegated, but via pre-signed URLs the Platform returns, so no cloud provider SDK integration is required.
- Q: How are credentials handled for the cloud object storage behind a data-link? → A: Platform-brokered. The user provides only the `tower.accessToken`. The Platform returns a short-lived pre-signed URL for each read; no cloud credentials cross the plugin boundary.
- Q: Read-only, or read + write for data-links in this iteration? → A: Read-only. Write (upload) support may be added later; the architecture does not preclude it but the feature is not in scope.
- Q: Path hierarchy — does the data-link identity segment include the provider? → A: Yes. `data-links/<provider>/<name>/...`. Names are not globally unique within a workspace (the same name may exist on two different providers), so the provider segment is required to disambiguate and mirrors the Platform UI's provider-grouped data explorer.
- Q: How deep can a data-link path go? → A: Arbitrary depth below the data-link root. Each segment after `<name>` is an entry inside the underlying bucket/prefix — a directory or file, resolved via the Platform browse API.
- Q: How should the existing dataset filesystem code be extended to accommodate data-links? → A: Introduce a true resource-type abstraction (`ResourceTypeHandler`). The current dataset-specific logic in `SeqeraFileSystemProvider`, `SeqeraFileSystem`, and `SeqeraPath` is extracted into a `DatasetsResourceHandler`; data-links are added as a parallel `DataLinksResourceHandler`. The core path/filesystem/provider classes become resource-type-agnostic.
- Q: How should the listing vs I/O boundary work? → A: Listing (`newDirectoryStream`) and attributes (`readAttributes`) are resolved via `GET /data-links/{id}/content?path=<sub>`. Downloads (`newInputStream`) go through `GET /data-links/{id}/download?path=<sub>`, which returns a pre-signed URL that is then fetched with the existing `TowerClient.sendStreamingRequest()` streaming helper. No cloud SDK is used.
- Q: Which DTOs are introduced by this feature? → A: None. All types are reused from the `io.seqera:tower-api:1.121.0` dependency (`DataLinkDto`, `DataLinkItem`, `DataLinkProvider`, `DataLinkContentResponse`, `DataLinkDownloadUrlResponse`, etc.).
- Q: Is browse-per-file supported by the Platform API? → A: Yes. `GET /data-links/{id}/content?path=<sub>` works for both directories and files, so `readAttributes` on any path is a single targeted call — no parent-browse-and-filter, no N+1 problem.
- Q: What happens if the pre-signed URL expires during a long read? → A: The underlying HTTP connection errors out with an `IOException`. The plugin does not transparently re-issue URLs; Nextflow's task retry handles the failure as it already does for other transient I/O errors.

## User Scenarios & Testing *(mandatory)*

### User Story 1 — Use a File Inside a Data-Link as Pipeline Input (Priority: P1)

A Nextflow pipeline developer has registered an S3 bucket or a GCS prefix as a Seqera Platform data-link (e.g. `inputs` on AWS). They want to reference a file inside that data-link as a pipeline input using a `seqera://` path, without configuring cloud credentials separately and without any manual pre-download step.

**Why this priority**: This is the core value proposition of the feature. All other stories build on the ability to resolve and read a file inside a data-link by path.

**Independent Test**: Write a pipeline that sets an input channel to `seqera://<org>/<ws>/data-links/<provider>/<name>/<relative-path>` and verify the pipeline task receives the correct file content, using only the Seqera access token for authentication.

**Acceptance Scenarios**:

1. **Given** a data-link named `inputs` registered with provider `aws` in workspace `acme/research`, pointing at `s3://my-bucket/data/`, **When** a pipeline references `seqera://acme/research/data-links/aws/inputs/reads/sample1.fq.gz`, **Then** the pipeline task receives the byte content of `s3://my-bucket/data/reads/sample1.fq.gz` transparently.
2. **Given** a data-link on `google` and one on `azure` in the same workspace, **When** both are referenced from the same pipeline, **Then** each path resolves independently and content is streamed correctly from the respective provider.
3. **Given** a Seqera access token is configured (`tower.accessToken` or `TOWER_ACCESS_TOKEN`), **When** a data-link path is accessed, **Then** no additional cloud credentials (AWS, GCP, Azure) are required from the user.

---

### User Story 2 — Browse the Data-Link Hierarchy (Priority: P2)

A pipeline developer wants to navigate the data-link namespace at any level — list providers in a workspace, list data-links within a provider, and browse into the content tree of a specific data-link — using ordinary directory listing operations.

**Why this priority**: Hierarchical listing supports discoverability and dynamic pipeline construction. Without it, users must know the exact path in advance.

**Independent Test**: List each level of the hierarchy and verify the correct child entries appear.

**Acceptance Scenarios**:

1. **Given** a workspace `acme/research` has data-links on both AWS and GCS, **When** a user lists `seqera://acme/research/data-links/`, **Then** `aws` and `google` are returned as directory entries (only providers in use appear).
2. **Given** the `aws` provider has two data-links `inputs` and `archive`, **When** a user lists `seqera://acme/research/data-links/aws/`, **Then** both names are returned as directory entries.
3. **Given** a data-link `inputs` contains a folder `reads/` with files `a.fq.gz` and `b.fq.gz`, **When** a user lists `seqera://acme/research/data-links/aws/inputs/reads/`, **Then** both files appear as file entries with correct size and last-modified metadata.
4. **Given** a data-link root is listed, **When** the data-link is empty, **Then** an empty result is returned without errors.
5. **Given** a user lacks access to a workspace, **When** they attempt to list any `data-links/` path within it, **Then** a clear access-denied error is returned without leaking internal details.
6. **Given** `readAttributes` is called directly on a file path (no prior listing), **When** the path exists, **Then** a single Platform API call returns the file attributes (no parent-directory scan).

---

### User Story 3 — Receive Meaningful Errors for Invalid or Inaccessible Paths (Priority: P3)

A pipeline developer mistypes a data-link name, uses a provider segment that has no registered data-links, or references a path that no longer exists inside a data-link. They receive a clear, actionable error that helps them fix the problem.

**Why this priority**: Error handling is essential for usability but delivers no new functionality on its own. Good errors prevent support escalations.

**Independent Test**: Reference invalid data-link paths and verify the error messages identify the problem (unknown provider, unknown data-link name, path not found inside data-link, authentication failure) without generic or cryptic failures.

**Acceptance Scenarios**:

1. **Given** a workspace has no data-links on the `azure` provider, **When** a pipeline references `seqera://.../data-links/azure/anything`, **Then** a `NoSuchFileException` is raised with a message indicating no data-links exist for provider `azure` in the workspace.
2. **Given** a data-link name that does not exist under a provider, **When** a pipeline attempts to read any path within it, **Then** the error identifies the missing data-link by name and path.
3. **Given** a valid data-link but a path that does not exist inside it, **When** a pipeline attempts to read it, **Then** a `NoSuchFileException` is raised with a message including the sub-path.
4. **Given** invalid or expired credentials, **When** any data-link path is accessed, **Then** an authentication error is reported with guidance to reconfigure `tower.accessToken` / `TOWER_ACCESS_TOKEN`.
5. **Given** a resource-type segment that is neither `datasets` nor `data-links`, **When** used, **Then** a clear "unsupported resource type" error is returned (unchanged from dataset feature).

---

### User Story 4 — Extensible Resource-Type Architecture (Priority: P4)

A Nextflow or Seqera engineer wants the filesystem's resource-type abstraction to be real and exercised by more than one resource type, so that adding future Seqera-managed resources requires isolated, scoped changes rather than cross-cutting refactors.

**Why this priority**: Shipping a second resource type is the right moment to introduce the shared abstraction — with two concrete consumers in place, the interface is validated in practice. This story captures the refactor as a first-class deliverable alongside the new data-link functionality.

**Independent Test**: A code review confirms that (a) `SeqeraFileSystemProvider`, `SeqeraFileSystem`, and `SeqeraPath` contain no dataset- or data-link-specific branching for depth ≥ 3; (b) adding a hypothetical third resource type is a new `ResourceTypeHandler` implementation with no changes to the generic path/provider/filesystem classes; (c) both `DatasetsResourceHandler` and `DataLinksResourceHandler` implement the same interface without leakage of each other's concepts.

**Acceptance Scenarios**:

1. **Given** the refactored filesystem, **When** `seqera://<org>/<ws>/` is listed, **Then** the resource-type entries (`datasets`, `data-links`) are enumerated from the handler registry rather than a hard-coded list.
2. **Given** the refactored `SeqeraPath`, **When** any path shape valid for either resource type is parsed, **Then** parsing succeeds without requiring the path class to know which resource type owns the depth-4+ segments.
3. **Given** a new handler is registered in `SeqeraFileSystem`, **When** paths with its resource-type segment are resolved, **Then** dispatch reaches it without modifying existing handlers.
4. **Given** both existing resource types, **When** a path from one is accessed, **Then** the other handler's code is never executed.

---

### Edge Cases

- What happens when a data-link's underlying bucket/prefix has been revoked on the provider side but the data-link still exists in Seqera? The Platform surfaces an error which is propagated as `IOException`.
- What happens when a data-link has thousands of entries at one level? The Platform browse endpoint's pagination (if any) must be exhausted; initial implementation pages through using whatever cursor token the response exposes.
- What happens when a user has a provider name containing characters the path class rejects? Provider identifiers come from the Platform API verbatim (`DataLinkProvider` enum values); they are valid path segments by construction.
- What happens when the same data-link name exists for two providers (e.g. `inputs` on AWS and `inputs` on GCS)? Both are addressable at distinct paths: `data-links/aws/inputs/...` and `data-links/google/inputs/...`. The provider segment disambiguates.
- What happens when a pre-signed URL expires mid-read? The HTTP read fails with `IOException`. The plugin does not transparently re-issue URLs; Nextflow task retry handles the failure.
- What happens when the transient Platform API is unavailable? Same as the dataset feature — `TowerClient`'s retry/backoff is reused; exhaustion raises `IOException`.
- What happens when a data-link's listing contains entries whose names include `/`? The browse response is expected to return name segments, not paths; any entry whose name contains `/` is rejected with a descriptive error (indicates a provider data issue).
- What happens when a data-link is accessed concurrently from many pipeline tasks? All reads are independent signed-URL fetches; no shared state beyond the (read-only) cached data-link list.
- How is pagination of `GET /data-links` handled if a workspace has more data-links than fit in a single response page? Implementation must exhaust pages before caching.

## Requirements *(mandatory)*

### Functional Requirements

- **FR-001**: System MUST accept paths in the format `seqera://<org>/<workspace>/data-links/<provider>/<name>/<sub-path>` where `<sub-path>` is zero or more segments addressing a directory or file inside the data-link.
- **FR-002**: System MUST read file content addressed by a data-link path transparently, requiring only the existing `tower.accessToken` / `TOWER_ACCESS_TOKEN` configuration — no cloud-provider credentials.
- **FR-003**: System MUST perform listing and attribute queries via the Seqera Platform API (`GET /data-links/{id}/content?path=<sub>`) and stream file content via pre-signed URLs returned from `GET /data-links/{id}/download?path=<sub>`.
- **FR-004**: System MUST support hierarchical directory listing:
  - `seqera://<org>/<ws>/` → directory; entries include `datasets` and `data-links` (enumerated from the handler registry).
  - `seqera://<org>/<ws>/data-links/` → directory; entries are distinct provider identifiers present in the workspace.
  - `seqera://<org>/<ws>/data-links/<provider>/` → directory; entries are data-link names under that provider.
  - `seqera://<org>/<ws>/data-links/<provider>/<name>/` → directory; entries are the top-level items in the data-link.
  - `seqera://<org>/<ws>/data-links/<provider>/<name>/<sub-path>/` → directory; entries are the children at that sub-path.
  - `seqera://<org>/<ws>/data-links/<provider>/<name>/<sub-path>/<file>` → file.
- **FR-005**: System MUST return correct `BasicFileAttributes` — `isDirectory`, `isRegularFile`, `size`, `lastModifiedTime`, `creationTime` — for any path inside a data-link, sourced from the Platform's content response for that path.
- **FR-006**: System MUST treat data-link paths as read-only in this iteration. Any write-like operation (`newByteChannel` with `WRITE`/`APPEND`, `copy` with a data-link as target, `delete`, `createDirectory`, `move`) MUST fail with `UnsupportedOperationException` or `AccessDeniedException`, consistent with the dataset feature's read-only stance.
- **FR-007**: System MUST produce clear, actionable error messages distinguishing: unknown org/workspace, unknown provider, unknown data-link name, missing sub-path, unsupported resource type, authentication failure, and transient Platform errors.
- **FR-008**: System MUST NOT depend on `nf-amazon`, `nf-google`, or `nf-azure`. All cloud I/O is reduced to a single HTTPS fetch of a pre-signed URL via the existing `TowerClient` / `HxClient` streaming path.
- **FR-009**: System MUST reuse DTOs from `io.seqera:tower-api:1.121.0` (`DataLinkDto`, `DataLinkContentResponse`, `DataLinkItem`, `DataLinkDownloadUrlResponse`, `DataLinkProvider`, etc.) without introducing parallel plugin-local classes.
- **FR-010**: System MUST refactor the existing `fs/` package to introduce a `ResourceTypeHandler` interface. `DatasetsResourceHandler` MUST encapsulate all dataset-specific behavior previously inlined in `SeqeraFileSystemProvider` / `SeqeraFileSystem` / `SeqeraPath`. `DataLinksResourceHandler` MUST implement the same interface.
- **FR-011**: After the refactor, the classes `SeqeraPath`, `SeqeraFileSystem`, and `SeqeraFileSystemProvider` MUST contain no dataset- or data-link-specific logic for depth ≥ 3; all such logic MUST live in the respective handler.
- **FR-012**: `SeqeraPath` MUST parse and represent arbitrary sub-paths below depth 4 for resource types that support them (data-links). Datasets continue to reject sub-paths beyond depth 4.
- **FR-013**: The filesystem MUST reuse the existing `TowerClient` retry/backoff for all Platform API calls. No new retry logic is introduced.
- **FR-014**: Transient failure of a pre-signed URL fetch mid-stream MUST surface as `IOException`; Nextflow task retry handles the recovery. The plugin MUST NOT re-issue URLs transparently within a single `InputStream`.
- **FR-015**: Data-link list results (per workspace) MUST be cached for the lifetime of a single `SeqeraFileSystem` instance. No browse-result or URL cache is maintained.
- **FR-016**: `GET /data-links?workspaceId=X` pagination MUST be exhausted before the workspace data-link list cache is considered populated.

### Key Entities

- **Data-Link**: A Seqera Platform entity referencing a bucket or prefix on a cloud provider (S3, GCS, Azure Blob, etc.). Addressed by `(workspaceId, provider, name)`; content is browsed and read through Platform API calls. Represented in the path as `data-links/<provider>/<name>/<sub-path>`.
- **Data-Link Provider**: A Platform-defined identifier for the cloud backend (`DataLinkProvider` enum values, e.g. `aws`, `google`, `azure`). Used as a path segment to disambiguate data-links with the same name on different providers.
- **Data-Link Entry**: An item inside a data-link — a file or folder — returned by the content API. Has a name, type (`FILE`/`FOLDER`), size, and last-modified timestamp.
- **Resource-Type Handler**: A pluggable strategy that owns the semantics of one depth-3 path segment (`datasets`, `data-links`, …). Exposes listing, attribute, read, and access-check operations to the generic filesystem.
- **Seqera Path (data-link variant)**: The URI `seqera://<org>/<ws>/data-links/<provider>/<name>[/<sub>/…]`. All segments up to and including `<name>` form the data-link identity; the remainder is the sub-path within the data-link.

## Success Criteria *(mandatory)*

### Measurable Outcomes

- **SC-001**: A pipeline developer can reference a file inside a Seqera data-link via `seqera://` path, and the pipeline runs successfully using only the Seqera access token — no cloud credentials, no manual pre-download step.
- **SC-002**: 100% of existing Nextflow file operations that work on cloud-hosted files (read, iterate lines, pass as channel input) work identically when the file is referenced via a `seqera://` data-link path.
- **SC-003**: Listing any level of the data-link hierarchy completes in under 5 seconds for workspaces with up to 500 data-links and data-links containing up to 1,000 entries at a single level.
- **SC-004**: Invalid or inaccessible data-link paths produce error messages that allow a developer to identify and fix the problem without consulting external documentation in 90% of cases (measured by user testing or code review of error-text coverage).
- **SC-005**: The refactored `fs/` package passes a code review confirming that (a) no resource-type-specific logic remains in `SeqeraPath`, `SeqeraFileSystem`, or `SeqeraFileSystemProvider`; (b) the two handlers share no code paths for depth ≥ 3; (c) the dataset tests pass unchanged after routing through `DatasetsResourceHandler`.
- **SC-006**: The plugin's runtime classpath for this feature gains no new cloud-SDK dependency (no `aws-sdk`, no `google-cloud-storage`, no `azure-*` artifacts introduced by this change).

## Assumptions

- Authentication reuses the existing nf-tower plugin credential mechanism (Seqera access token); no new auth configuration is required from users.
- The `GET /data-links/{id}/content?path=<sub>` endpoint works for both directory and file paths. When `path` points at a file, the response describes just that file; when it points at a directory, the response's `items` array enumerates children.
- The `GET /data-links/{id}/download?path=<sub>` endpoint returns a pre-signed URL valid for long enough to complete a typical file read. The plugin does not extend this window.
- Data-link provider identifiers returned by the Platform (`DataLinkProvider`) are safe as path segments — they contain no `/` and no characters that collide with URI reserved characters.
- The tower-api artifact (`io.seqera:tower-api:1.121.0`) already available on the plugin classpath exposes all DTOs required (`DataLinkDto`, `DataLinkContentResponse`, `DataLinkItem`, `DataLinkDownloadUrlResponse`, `DataLinkProvider`, etc.).
- Data-link writes, renames, deletes, and management operations (create, update, delete the data-link entity itself) are **out of scope** for this iteration.
- Streaming a pre-signed URL reuses `TowerClient.sendStreamingRequest()` and therefore inherits the same retry/backoff behavior as dataset downloads.
- Data-link listings may be paginated; the plugin exhausts pagination before caching.
- No local caching across pipeline runs. Nextflow's standard task staging handles intra-run caching.
- Paths are case-sensitive — matches the Platform API and the dataset filesystem.
- The dataset feature's read-only filesystem stance (`isReadOnly()=true`) is preserved; data-link writes are deferred to a future iteration.

## Dependencies

- Seqera platform API (data-links endpoints: `/data-links`, `/data-links/{id}/content`, `/data-links/{id}/download`) must be accessible from the compute environment where the pipeline runs.
- nf-tower plugin must be enabled and configured with a valid `tower.accessToken` / `TOWER_ACCESS_TOKEN`.
- The Seqera account must have at least read access to the target workspace and data-link.
- The existing dataset filesystem (`260310-seqera-dataset-fs`) must be merged — this feature builds on its classes and refactors them.

## Out of Scope

- Write operations (upload) to data-links — the Platform's `POST /data-links/{id}/multipart-upload` endpoint is a natural future hook but is not implemented here.
- Data-link management operations (create, update, delete the data-link entity itself).
- Transparent pre-signed URL renewal mid-stream.
- Local caching across pipeline runs.
- Browse-result caching within a run.
- Fusion integration (Fusion has its own data-link access path; this feature is for direct NIO access).
