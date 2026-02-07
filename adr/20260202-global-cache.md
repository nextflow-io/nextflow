# Global Cache for Cross-Pipeline Task Sharing

- Authors: Jorge Ejarque
- Status: draft
- Deciders: Paolo Di Tommaso, Ben Sherman, Phil Ewels
- Date: 2026-02-02
- Tags: cache, performance, task-reuse, cross-pipeline, cloud-storage, nf-cloudcache
- Version: 1.0 

## Summary

Extend Nextflow's task caching to enable cross-pipeline task result sharing through a global cache backed by cloud object storage (S3, GCS, Azure Blob), building on the existing `nf-cloudcache` infrastructure and leveraging cloud storage strong consistency guarantees for concurrent access control.

## Problem Statement

Nextflow's current caching mechanism (`-resume`) is designed exclusively for resuming the same pipeline execution after interruption. This design creates significant limitations:

1. **Pipeline-specific hashing**: Task hashes include the process name and session ID, making them unique to each pipeline execution. Different pipelines computing identical tasks cannot reuse each other's results.

2. **File hashing limitations**: Depending on the cache mode, file hashing may include the file path and filesystem attributes instead of actual content, preventing cache hits when the same data exists at different paths.

3. **No cross-pipeline sharing**: Organizations running similar or related pipelines repeatedly compute the same intermediate results, wasting compute resources and time.

4. **Session-bound concurrency control**: The current `LockManager` (TaskProcessor.groovy:825-839) prevents race conditions within a single session but doesn't coordinate across different pipeline executions or machines.

These limitations result in redundant computation when:
- Different users run the same analysis on the same datasets (different caches)
- Similar pipelines share common processing steps
- Development and production versions of a pipeline process identical data
- Parametric sweeps or batch analyses repeat computations
- Multiple pipeline executions occur across different machines/environments

## Goals or Decision Drivers

- **Cross-pipeline task reuse**: Enable different pipelines to reuse task results when inputs and processing logic are identical
- **Content-addressable hashing**: Hash tasks based on actual content rather than pipeline-specific metadata (no process name, no session ID)
- **Cloud-native backend**: Leverage existing `nf-cloudcache` infrastructure and cloud storage (S3, GCS, Azure Blob)
- **Strong consistency guarantees**: Use cloud storage consistency for concurrent access control across pipeline executions
- **Simplicity and reliability**: Prefer simple collision-avoidance over complex waiting/polling mechanisms (following existing Nextflow patterns)
- **Minimal new infrastructure**: Build on existing cloud storage used for work directories
- **Backward compatibility**: Maintain existing cache behavior as default
- **Transparent adoption**: Minimal changes to existing pipelines

## Non-goals

- **Maintaining local filesystem cache**: Global cache is cloud storage only
- **Complex waiting/polling**: No waiting for in-progress tasks; use hash increment instead (simplicity over optimal cache hits)
- **Perfect concurrent coordination**: Accept rare redundant execution during race conditions
- **Cache versioning**: No complex versioning or lineage tracking
- **Semantic caching**: No attempt to understand task equivalence beyond hash matching
- **Shared filesystem support**: Focus on cloud object storage, not NFS/Lustre
- **Automatic migration**: No automatic migration of existing local cache

## Solution Approach

Extend the existing `nf-cloudcache` plugin to support content-addressable global caching on cloud object storage (S3, GCS, Azure Blob).

**Rationale for this approach:**
- `nf-cloudcache` already exists and handles cloud storage integration
- Cloud storage provides strong consistency guarantees for concurrent access
- Many organizations already use cloud storage for work directories
- Cloud providers support atomic operations needed for coordination
- No new infrastructure required
- Scalable and accessible from anywhere

**Trade-offs:**
- Higher latency than local filesystem (~100-500ms vs ~10ms)
- Cloud storage API costs (offset by compute savings from cache hits)
- Requires cloud credentials (already needed for cloud executors)

## Alternative Approaches (For Reference)

Alternative approaches that could be considered in the future but are not part of this design:

**Shared Filesystem (NFS/Lustre):**
- Not pursued: Requires dedicated infrastructure, weaker consistency, limited to single site

**Hybrid Local + Global Cache:**
- Not pursued: Adds complexity with two cache systems and coherence issues

**Distributed Cache Service (Redis, etc.):**
- Not pursued: Requires new infrastructure, doesn't scale for large outputs, operational burden

## Rationale & discussion

### Architecture Overview

Build on existing `nf-cloudcache` infrastructure:

```
┌─────────────────────────────┐
│   Task Execution Engine     │
│   (TaskProcessor.groovy)    │
└──────────────┬──────────────┘
               │
      Check global cache
               │
    ┌──────────▼──────────┐
    │  Global Cache       │
    │  (nf-cloudcache     │
    │   extension)        │
    └──────────┬──────────┘
               │
       ┌───────┴────────┐
       │                │
   Cache hit        Cache miss
       │                │
   Reuse result    Execute task
       │                │
       └────────┬───────┘
                │
         Store in cloud
         (S3/GCS/Azure)
```

**Cache lookup flow:**
1. Compute global hash based on task content (no sessionId, no processName)
2. Check cloud cache for matching hash prefix
3. If workdir and `.exitcode` exists and valid: cache hit, restore outputs
4. If workdir exists and `.exitcode` not being created: task in progress. Increment hash and try again
5. If no match: cache miss, execute task normally

### Hash Algorithm Redesign

**Current local cache hash (TaskProcessor.groovy):**
```groovy
localHash = hash(
    sessionId,           // ❌ Prevents cross-session sharing
    processName,         // ❌ Prevents cross-pipeline sharing
    script,
    inputs,              // Includes file paths/attributes
    container,
    environment
)
```

**New global cache hash:**
```groovy
globalHash = hash(
    script,              // Task script content
    inputsContent,       // ✓ Content-based input hashing (checksums)
    containerDigest,     // ✓ Container image SHA256 digest (not tag)
    environment
)
```

**Key differences:**
- **No sessionId**: Enables sharing across all sessions and pipeline executions
- **No processName**: Enables sharing across different process names with identical logic
- **Content-based inputs**: Files hashed by full content (checksum), not path/metadata
- **Container digest**: Use SHA256 digest (e.g., `sha256:abc123...`) instead of tag

### Input File Hashing Alternatives

The global cache requires a method to fingerprint input files as part of the task hash computation. Unlike local resume cache where file paths remain stable across executions, global cache must handle scenarios where identical file content exists at different paths, locations, or across different pipeline executions. The choice of hashing strategy directly impacts cache hit rates, performance overhead, and cross-execution compatibility.

#### Alternative 1: File Path and Attributes (Current Default for Local Resume)

**Description:** Hash based on file path, size, and modification timestamp.

**Pros:**
- Very fast (~1ms per file)
- Minimal computational overhead
- Works perfectly for local resume with stable paths

**Cons:**
- Different paths produce different hashes even for identical content
- Timestamps and paths change across environments

**Cache hit scenarios in global cache:**
- **Low hit rate** - only when same absolute paths used:
  - Shared filesystem with consistent mount points (e.g., `/data/refs/hg38.fa` always at same path)
  - Published outputs accessed from fixed publish directories
  - Standardized organizational data layouts
- **Misses when:**
  - Same file copied/moved to different path
  - Files accessed from different machines/mount points
  - Cross-user sharing (different home directories)

**Performance:** ⚡ Fastest - negligible overhead

---

#### Alternative 2: Deep Content Hashing (SHA256 of Full File)

**Description:** Compute cryptographic hash of entire file content.

**Pros:**
- Content-based - identical content always produces same hash regardless of path
- Cryptographically secure and collision-resistant
- Enables true cross-pipeline and cross-environment sharing
- No external dependencies required

**Cons:**
- High computational cost (~100-500 MB/s throughput)
- Significant overhead for large files (10 GB file = ~20-100 seconds)
- Repeated computation for same file used by multiple tasks
- I/O intensive on network storage

**Cache hit scenarios in global cache:**
- **High hit rate** - whenever file content matches:
  - Same content at different paths (e.g., `/user1/genome.fa` vs `/user2/refs/hg38.fa`)
  - Files copied, moved, or republished
  - Cross-pipeline intermediate result sharing
  - Reference data accessed from various locations
- **Always achieves hit for identical content** (correctness guaranteed)

**Best for:**
- Small to medium files (<100 MB)
- High-value computations where cache hit saves significant time
- Scenarios where correctness is critical

**Performance:** Slowest - O(file size)

---

#### Alternative 3: Lineage Metadata-Based (`lid` Fingerprinting)

**Description:** Use Nextflow lineage tracking to identify file content via lineage ID (`lid`). Check if lineage entry exists for file path; if yes, use original `lid` as fingerprint; if no, fall back to deep hash and create lineage entry.

**Pros:**
- Fast lookup for tracked files (~1ms)
- Avoids re-hashing Nextflow-managed files
- Maintains content-based semantics
- Leverages existing lineage infrastructure

**Cons:**
- Requires extending lineage schema to track content hashes
- Lineage database must be shared across executions
- Only works for Nextflow-managed files
- External files require fallback to deep hash

**Cache hit scenarios in global cache:**
- **High hit rate for Nextflow ecosystem files:**
  - Task output reused by downstream tasks (lineage tracks creation)
  - Published files consumed by other pipelines (lineage tracks publish)
  - Same file referenced from multiple paths but same `lid`
- **Requires deep hash fallback for:**
  - External reference data not created by Nextflow (first use only)
  - User-uploaded files (first use only)
  - Files from non-Nextflow sources
- **Example:** Preprocessing output used by 50 downstream tasks - first task deep hashes and creates `lid`, remaining 49 use fast lookup

**Best for:**
- Workflows with many intermediate file reuse patterns
- DAG fan-out (one file → many consumers)
- Environments where most inputs are Nextflow-generated outputs

**Limitations:**
- Only achieves cache hits for files previously used or published by Nextflow within the global cache system
- Files from external sources (uploaded directly to cloud, other tools) cannot be matched even if content identical
- Most effective in closed ecosystem where Nextflow manages all files

**Performance:** Variable: Fast for tracked files, Slow for deep hash fallback for external files

**Implementation requirements:**
- Extend lineage schema to store content hash alongside `lid`
- File path → lineage entry index for lookups
- Lineage validation (verify file unmodified since tracking)
- Lineage database synchronization across distributed executions

---

#### Alternative 4: Cloud Storage Provider Checksums

**Description:** Use checksums computed by cloud storage (S3 ETags, GCS MD5, Azure MD5) instead of local computation.

**Pros:**
- Fast (~100ms metadata API call vs seconds for deep hash)
- No local computation - provider pre-computes checksums
- Works for any file in cloud storage

**Cons:**
- **No consistent algorithm across providers:**
  - GCS/Azure: MD5 (not collision-resistant)
  - S3: Multiple algorithms depending on upload method
- **MD5 security concerns:** Vulnerable to collisions
- **S3 inconsistency:** Different objects may use MD5, SHA256, CRC32
- **Multipart uploads:** S3 ETags are composite hashes, not standard checksums

**Cache hit scenarios in global cache:**
- **Medium hit rate with significant constraints:**
  - **GCS or Azure only deployments:** Consistent MD5 across all objects → high hit rate within provider, but MD5 collision risk
  - **S3 with consistent upload method:** Hits when same algorithm used (e.g., all files uploaded with SHA256 via Nextflow)
- **Misses when:**
  - Same content uploaded with different algorithms in S3 (one upload MD5, another SHA256)
  - Files transferred between cloud providers (GCS MD5 ≠ S3 SHA256)
  - Multipart vs single-part uploads in S3 (different ETag formats)
  - External files uploaded without checksums or with incompatible algorithms
- **Example hits:** Pure Azure deployment where all files use MD5, or S3 environment with governance policy enforcing SHA256 for all uploads
- **Example misses:** File uploaded to S3 as multipart (composite ETag) cannot match same content uploaded as single object (standard checksum)

**Provider-specific characteristics:**

| Provider | Algorithm | Collision-Resistant | Consistent | Multipart Handling | Reliability |
|----------|-----------|--------------|------|-------------------|-------------|
| GCS | MD5 | No | Yes | Standard MD5 | Medium (consistent but MD5 risk) |
| Azure | MD5 | No | Yes | Standard MD5 | Medium (consistent but MD5 risk) |
| S3 | Varies | Depends | No | Composite ETag | Low (algorithm inconsistency) |

**Detailed provider constraints:**

- **Google Cloud Storage:** Always MD5. Consistent across all objects, enabling reliable matching within GCS. However, MD5 is not cryptographically secure and vulnerable to collisions.

- **Azure Blob Storage:** Always MD5. Similar to GCS - consistent within Azure, but MD5 collision vulnerability.

- **AWS S3:** Supports multiple algorithms (MD5, SHA256, CRC32, CRC32C). Algorithm varies based on:
  - Upload method (PUT vs multipart)
  - Client SDK settings and version
  - PUT request headers (can specify SHA256 via `x-amz-checksum-sha256`, see #6321)
  - Historical objects uploaded before SHA256 support
  - **Multipart uploads use composite ETags** (hash of chunk hashes, not file content hash)

  **Result:** Cannot reliably match files across different upload methods or time periods within same S3 bucket.

**Best for:**
- Proof-of-concept implementations
- Trusted environments without collision attack concerns
- Single cloud provider with enforced upload standards
- Performance-critical scenarios accepting correctness trade-offs

**Performance:** Fast - ~100ms API call

---

#### Alternative 5: Hybrid Tiered Approach

**Description:** Combine multiple methods with fallback chain:
1. Lineage metadata (if available) → fast
2. Cloud provider checksums (if single-provider environment) → medium
3. Deep content hash (for external/mixed data) → slow but correct

**Pros:**
- Optimizes common cases (lineage) while ensuring correctness (deep hash fallback)
- Adaptable to different deployment scenarios
- Balances performance and correctness

**Cons:**
- Implementation complexity with multiple code paths
- Difficult to reason about behavior and cache hit rates
- Requires lineage database synchronization
- Inconsistent performance depending on file source

**Cache hit scenarios:**
- Variable depending on which tier matches
- Best case: Lineage hit (fast, high reliability)
- Medium case: Cloud checksum hit (medium speed, medium reliability)
- Worst case: Deep hash (slow, high reliability)

**Best for:**
- Production deployments with diverse file sources
- Mixed workloads (internal + external data)
- Organizations willing to accept complexity for optimization

**Performance:** Variable - 1ms to minutes depending on tier

---

#### Comparison Summary

| Alternative | Performance | Cache Hit Rate (Global) | Correctness | Implementation Complexity | Best Environment |
|-------------|-----------|------------------------|-------------|--------------------------|------------------|
| Path/Attributes | ⚡⚡⚡ Fastest | 🔴 Low | ❌ Path-dependent | ✅ Simple | Local resume only |
| Deep Content Hash | Slowest | 🟢 Highest | ✅ Guaranteed | ✅ Simple | Universal (MVP choice) |
| Lineage Metadata | ⚡⚡ Fast | 🟢 High (Nextflow files) | ✅ Good | ⚠️ Medium | Closed Nextflow ecosystem |
| Cloud Checksums | ⚡⚡ Fast | 🟡 Medium | ⚠️ Provider-dependent | ⚠️ Medium | Single-provider, trusted |
| Hybrid Tiered |  Variable | 🟢 High (variable) | ✅ Good | 🔴 Complex | Production, mixed workloads |



### Concurrent Access Control

**The Problem:**

The current mechanism (TaskProcessor.groovy:825-839) checks if work directory exists:
```groovy
final lock = lockManager.acquire(hash)  // Only coordinates within one session
try {
    if( workDir.exists() ) {
        tries++
        continue  // Directory exists, try next hash
    }
    else if( !workDir.mkdirs() ) {
        throw new IOException("Unable to create directory")
    }
}
finally {
    lock.release()
}
```

**Issue with global cache:**
- `LockManager` only coordinates within a single Nextflow session
- Two pipelines on different machines can both check the same cache directory doesn't exist
- Both pipelines try to create it simultaneously → race condition
- Without atomic operations, both could create the directory and execute redundantly
- Strict consistency do not solve the issue by itself because two put operation could happen at the same time. (require preconditions to produce a failure it it happens)

**Solution:** Use cloud storage atomic operations with preconditions.

#### Cloud Provider Atomic Operations

**AWS S3 - Conditional PUT:**
```
PUT s3://bucket/cache/ab/cdef123.../.lock
Headers:
  If-None-Match: *    # Only create if object doesn't exist

Response:
  200 OK              # Success - we own this task
  412 Precondition Failed  # Object exists - another pipeline owns it
```

**Google Cloud Storage - Generation Preconditions:**
```
PUT gs://bucket/cache/ab/cdef123.../.lock
Parameters:
  ifGenerationMatch=0    # Generation 0 = object doesn't exist

Response:
  200 OK                 # Success - we own this task
  412 Precondition Failed  # Object exists - another pipeline owns it
```

**Azure Blob Storage - Conditional Headers:**
```
PUT https://account.blob.core.windows.net/cache/ab/cdef123.../.lock
Headers:
  If-None-Match: *    # Only create if blob doesn't exist

Response:
  201 Created          # Success - we own this task
  409 Conflict         # Blob exists - another pipeline owns it
```

**Note:** Azure also supports explicit blob leases for locking, but conditional headers are simpler and sufficient.

#### Strategy: Simple Collision Avoidance

Following Nextflow's existing philosophy (TaskProcessor.groovy:825-839), use **collision avoidance** instead of waiting:

**If atomic create fails (directory exists), increment hash and try next** (trade cache hit for simplicity).

**Cache lookup and execution flow:**
```
1. Compute globalHash from task content

2. tries = 0
3. Loop :
   a. currentHash = (tries == 0) ? globalHash : hash(globalHash, tries)

   b. Check if hash entry exists, <currentHash>/.exitcode exists:
      → YES: Cache hit! Restore outputs, done

   c. Try atomic create of <currentHash>/.lock with precondition:
      → SUCCESS: We own this cache entry, execute task
      → FAILURE: Someone else owns it, increment tries (

```

**Task execution after successful atomic create:**
```
5. Execute task in local work directory
6. Upload outputs to <currentHash>/
7. Upload .exitcode and add hash entry (marks completion)
```

**Why this works:**
- Atomic create prevents two pipelines from claiming same cache entry
- If directory exists (incomplete or in-progress), increment hash and try next
- Simple, no waiting/polling/timeouts
- Trade-off: Rare concurrent execution (~1% of tasks) executes redundantly

**Key Race Condition Scenarios:**

**Scenario 1: Simultaneous execution (most common)**
```
Time  Pipeline A                    Pipeline B
----  -------------------------     -------------------------
T0    Check hash123/                Check hash123/
      → missing                      → missing

T1    Atomic PUT hash123/.lock     Atomic PUT hash123/.lock
      → SUCCESS (200 OK)             → FAIL (412 Precondition Failed)

T2    Execute task                  Increment to hash124
                                    Try hash124 instead

Result: Both pipelines execute (redundant), but no data corruption
```

**Scenario 2: Crashed task left incomplete entry**
```
Time  Pipeline A (yesterday)        Pipeline B (today)
----  -------------------------     -------------------------
T0    Atomic PUT hash123/.begin
      → SUCCESS

T1    Execute task
      CRASH (no .exitcode uploaded)

      [incomplete entry remains]

T2                                  Check hash123/.exitcode
                                    → missing

T3                                  Increment to hash124
                                    Try hash124 instead

Result: Pipeline B avoids incomplete entry
Note: Cleanup process removes old incomplete entries periodically
```

**Summary - Cloud-Specific Implementation:**

| Cloud Provider | Atomic Operation | Precondition Header/Param | Success Code | Conflict Code |
|---|---|---|---|---|
| **AWS S3** | Conditional PUT | `If-None-Match: *` | 200 OK | 412 Precondition Failed |
| **GCS** | Generation precondition | `ifGenerationMatch=0` | 200 OK | 412 Precondition Failed |
| **Azure Blob** | Conditional header | `If-None-Match: *` | 201 Created | 409 Conflict |


### Migration and Backward Compatibility

**Default behavior (no change):**
```groovy
// Existing pipelines continue using local cache
nextflow run pipeline.nf -resume
```

**Opt-in to global cache:**

```bash
# Set nf-cloudcache, workdir as subfolders of the global cache path. Activate resume and task hash with content.
export NXF_GLOCALCACHE_PATH=s3://bucket/global-cache 
nextflow run pipeline.nf
```

**No automatic migration:**
- Existing local cache entries are not migrated
- Global cache starts empty
- Over time, populates as tasks execute
- Old local cache can be cleaned up manually

### Consideration with other planned features: Automatic Workflow Cleanup

The global cache feature collides with the Automatic Workflow Cleanup feature. In that planned feature, intermediate files could be automatically removed if they are not required to compute the pipeline outputs or for resume anymore.
In a global cache, these intermediate files could be required in case of a cache hit for the tasks that produced these files.

### Performance Considerations

**Content hashing overhead:**
- Reading entire input files for checksums
- Possibly mitigated with proposed optimizations: Lineage and cloud storage checksums

**Cloud storage latency:**
- Higher latency than local disk (100-500ms vs <10ms)

**Cost considerations:**
- S3/GCS API request costs
- Storage costs for cached outputs
- **Benefit:** Reduced compute costs from cache hits typically far exceed storage costs


### Storage Management

**Cache growth:**
- Unbounded growth as more unique tasks are cached
- Monitor storage usage

**Incomplete entry cleanup:**
- Tasks that crash leave incomplete entries in the workdir (`.command.begin` present but no `.exitcode`)
- These entries prevent cache hits and waste storage
- Periodic cleanup required to remove stale incomplete entries

**Manual cleanup:**
```bash
# Clean all cache entries
nextflow cache clean -global
```

### Security and Access Control

**Trust model:**
- Global cache uses the bucket security and access control
- Any user with access to the bucket can read/write cache
- Cached outputs are visible to all who has access to the bucket



### Testing Strategy

**Unit tests:**
- Global hash computation (no sessionId, no processName)
- Content-based file hashing
- Hash increment logic (collision resolution)
- Atomic work directory creation

**Integration tests:**
- Multiple pipelines executing concurrently
- Race condition scenarios (simultaneous work directory creation)
- Cache hit/miss verification
- Incomplete entry handling

**Cloud provider tests:**
- S3 conditional PUT operations (if-none-match)
- GCS generation preconditions
- Azure blob conditional headers
- Strong consistency verification

**Performance benchmarks:**
- Hash computation overhead
- Cloud storage latency
- Cache hit rate improvements (with and without concurrency)
- Redundant execution frequency during race conditions
- Cost analysis

### Example Use Cases

**Use Case 1: Shared reference processing**
```groovy
process index_reference {
    input:
    path reference  // Same hg38.fa for all users but in different locations

    output:
    path "*.idx"

    script:
    """
    build_index ${reference}  
    """
}

// First user: Creates index, uploads to global cache
// Subsequent users: Instant cache hit, reuse index
```

**Use Case 2: Parameter sweeps**
```groovy
// Preprocess step is identical for all parameter values
process preprocess {
    input:
    path data

    output:
    path "preprocessed.dat"

    script:
    """
    expensive_preprocessing ${data}  // 1 hour
    """
}

params.values = [1, 2, 3, 4, 5, ..., 100]

// First run: preprocess executed once, cached
// Runs 2-100: instant cache hit
```

**Use Case 3: Dev/prod parity**
```groovy
// Development pipeline
process dev_analysis { ... }

// Production pipeline (identical logic, different name)
process prod_analysis { ... }

// Production reuses dev results (same script, same inputs)
```

### Known Limitations

**Initial implementation:**
- Cloud storage only (no shared filesystem)
- Manual cache lifecycle management
- Content hashing overhead for large files
- Storage costs for cached outputs
- **Concurrent execution loses cache hits**: If two pipelines simultaneously execute identical tasks, both will execute redundantly in different hashes and workdirs (trade-off for simplicity)

**Edge cases:**
- Non-deterministic tasks may cause issues (random number generation, timestamps in outputs)
- Container tags (not digests) may cause false cache hits with different images
- Very large outputs may exceed practical storage limits or cause large storage cost
- Incomplete cache entries from crashed tasks require periodic cleanup

### Implementation Plan

**Phase 0: Proof of concept (#6100)**
1. Associate nf-cloudcache path and workdir with the global-cache path and active resume by default
2. Constant sessionId (0000-000-000), remove processName from task hash
3. Optional: Use deep cache mode

**Phase 1: Core functionality*
1. Implement global hash algorithm (no sessionId, no processName)
2. Implement content-based file hashing

**Phase 2: Concurrency control**
1. Add cloud storage lock acquisition (S3 conditional PUT)
2. Test race condition handling

**Phase 4: Polish**
1. Add configuration options
2. Implement cache cleanup commands 
3. Documentation and examples

## Links

- [nf-cloudcache plugin](../plugins/nf-cloudcache/) - Foundation for global cache
- [CloudCacheConfig](../plugins/nf-cloudcache/src/main/nextflow/cache/CloudCacheConfig.groovy) - Configuration class
- [TaskProcessor.groovy](../modules/nextflow/src/main/groovy/nextflow/processor/TaskProcessor.groovy) - Cache checking logic (lines 825-839, 925-1001)
- [GridTaskHandler.groovy](../modules/nextflow/src/main/groovy/nextflow/executor/GridTaskHandler.groovy) - Exit code checking (lines 313-417)
- [S3 Strong Consistency](https://aws.amazon.com/s3/consistency/) - S3 consistency guarantees
- [GCS Consistency](https://cloud.google.com/storage/docs/consistency) - GCS consistency model
- [Wave containers](https://seqera.io/wave/) - Container provenance with SHA256 digests