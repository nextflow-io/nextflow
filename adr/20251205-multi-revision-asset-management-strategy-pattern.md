# Multi-Revision Asset Management with Strategy Pattern

- Authors: Jorge Ejarque
- Status: Approved
- Deciders: Jorge Ejarque, Ben Sherman, Paolo Di Tommaso
- Date: 2025-12-05
- Tags: scm, asset-management, multi-revision

## Summary

Nextflow's asset management system has been refactored to support multiple revisions of the same pipeline concurrently through a bare repository approach with shared object storage, while maintaining backward compatibility with legacy direct-clone repositories using the Strategy design pattern.

## Problem Statement

The original asset management system (`AssetManager`) cloned each pipeline directly to `~/.nextflow/assets/<org>/<project>/.git`, creating several limitations:

1. **No concurrent multi-revision support**: Only one revision of a pipeline could be checked out at a time, preventing concurrent execution of different versions
2. **Update conflicts**: Pulling updates while a pipeline was running could cause conflicts or corruption
3. **Testing limitations**: Users couldn't easily test different versions of a pipeline side-by-side

The goal was to enable running multiple revisions of the same pipeline concurrently (e.g., production on v1.0, testing on v2.0-dev) while maintaining efficient disk usage through object sharing.

## Goals or Decision Drivers

- **Concurrent multi-revision execution**: Must support running different revisions of the same pipeline simultaneously
- **Efficient disk usage**: Share Git objects between revisions to minimize storage overhead
- **Backward compatibility**: Must not break existing pipelines using the legacy direct-clone approach
- **API stability**: Maintain the existing `AssetManager` API for external consumers (K8s plugin, CLI commands, etc.)
- **Minimal migration impact**: Existing repositories should continue to work without user intervention
- **JGit compatibility**: Solution must work within JGit's capabilities to avoid relying on Git client installations
- **Atomic updates**: Downloading new revisions should not interfere with running pipelines

## Non-goals

- **Migration of existing legacy repositories**: Legacy repos continue to work as-is; no forced migration
- **Native Git worktree support**: Due to JGit limitations, not using Git's worktree feature
- **Revision garbage collection**: No automatic cleanup of old revisions (users can manually drop)
- **Multi-hub support**: Still tied to a single repository provider per pipeline

## Considered Options

### Option 1: Bare Repository with Git Worktrees

Use Git's worktree feature to create multiple working directories from a single bare repository.

**Implementation**:
- One bare repository at `~/.nextflow/assets/<org>/<project>/.git`
- Multiple worktrees at `~/.nextflow/assets/<org>/<project>/<revision>/`

- Good, because it's the native Git solution for multiple checkouts
- Good, because worktrees are space-efficient
- Good, because Git handles all the complexity
- **Bad, because JGit doesn't support worktrees** (deal-breaker)
- Bad, because requires native Git installation

**Decision**: Rejected due to JGit incompatibility 

### Option 2: Bare Repository + Clones per Commit + Revision Map File

Use a bare repository for storage and create clones for each commit, tracking them in a separate file.

**Implementation**:
- Bare repository at `~/.nextflow/assets/<org>/<project>/.nextflow/bare_repo/`
- Clones at `~/.nextflow/assets/<org>/<project>/.nextflow/commits/<commit-sha>/`
- Revision map file at `~/.nextflow/assets/<org>/<project>/.nextflow/revisions.json` mapping revision names to commit SHAs

- Good, because it works with JGit
- Good, because bare repo reduces remote repository interactions to checkout commits
- Good, because explicit revision tracking
- Bad, because disk space as git objects replicated in clones 
- Bad, because revision map file can become stale
- Bad, because requires file I/O for every revision lookup
- Bad, because potential race conditions on map file updates
- Bad, because adds complexity of maintaining external state

**Decision**: Initially implemented but later refined

### Option 3: Bare Repository + Shared Clones with Strategy Pattern

Similar to Option 2 but eliminate the separate revision map file by using the bare repository itself as the source of truth. Additionally, use the Strategy pattern to maintain backward compatibility with existing legacy repositories without requiring migration.

**Implementation**:
- Bare repository at `~/.nextflow/assets/.repos/<org>/<project>/bare/`
- Shared clones at `~/.nextflow/assets/.repos/<org>/<project>/commits/<commit-sha>/`
- Use bare repository refs to resolve revisions to commit SHAs dynamically
- JGit alternates mechanism for object sharing
- `AssetManager` as facade with unchanged public API
- `RepositoryStrategy` interface defining repository operations
- `LegacyRepositoryStrategy` for existing direct-clone behavior
- `MultiRevisionRepositoryStrategy` for new bare-repo approach
- Strategy selection based on environment variable or repository state detection

- Good, because no external state file to maintain
- Good, because bare repository is always in sync (fetched on updates)
- Good, because simpler and more reliable
- Good, because atomic updates (Git operations are atomic)
- Good, because works entirely within JGit
- Good, because zero migration needed for existing repositories
- Good, because maintains API compatibility
- Good, because allows gradual adoption
- Good, because isolates legacy code
- Good, because makes future strategies easy to add
- Neutral, because adds abstraction layer
- Bad, because requires resolution on every access (minimal overhead)
- Bad, because increases codebase size initially

**Decision**: Selected

## Solution or decision outcome

Implemented **Option 3 (Bare Repository + Shared Clones with Strategy Pattern)** for multi-revision support with backward compatibility. Multi-revision is the default for new repositories, while legacy mode is available via `NXF_SCM_LEGACY` environment variable.

## Rationale & discussion

### Multi-Revision Implementation

The bare repository approach provides efficient multi-revision support:

```
~/.nextflow/assets/.repos/nextflow-io/hello/
├── bare/                       # Bare repository (shared objects)
│   ├── objects/                # All Git objects stored here
│   ├── refs/
│   │   ├── heads/
│   │   └── tags/
│   └── config
│
└── commits/                    # Commit-specific clones
    ├── abc123.../              # Clone for commit abc123
    │   └── .git/
    │       ├── objects/        # (uses alternates → bare/objects)
    │       └── info/
    │           └── alternates  # Points to bare/objects
    │
    └── def456.../              # Clone for commit def456
        └── .git/

~/.nextflow/assets/nextflow-io/hello/
└── .git/                       # Legacy repo location (HYBRID state)
```

**Key mechanisms:**

1. **Bare repository as source of truth**: The bare repo is fetched/updated from the remote, keeping refs current
2. **Dynamic resolution**: Revisions (branch/tag names) are resolved to commit SHAs using the bare repo's refs
3. **Object sharing**: Clones use Git alternates to reference the bare repo's objects, avoiding duplication
4. **Atomic operations**: Each clone is independent; downloading a new revision doesn't affect existing ones
5. **Lazy creation**: Clones are created on-demand when a specific revision is requested

**Advantages over revision map file:**
- No external state to maintain or keep in sync
- Bare repo fetch automatically updates all refs
- Resolution is simple: `bareRepo.resolve(revision)` returns commit SHA
- No race conditions on file updates
- Simpler code with fewer failure modes

### Strategy Pattern for Backward Compatibility

The Strategy pattern provides clean separation and backward compatibility:

```
┌─────────────────────────┐
│     AssetManager        │  ← Public API (unchanged)
│     (Facade)            │
└───────────┬─────────────┘
            │
            │ delegates to
            ▼
┌─────────────────────────┐
│  RepositoryStrategy     │  ← Interface
└───────────┬─────────────┘
            △
            │ implements
    ┌───────┴────────┐
    │                │
┌───────────┐  ┌─────────────────┐
│  Legacy   │  │ MultiRevision   │  ← Concrete strategies
│ Strategy  │  │   Strategy      │
└───────────┘  └─────────────────┘
```

**Strategy selection logic:**

1. Check `NXF_SCM_LEGACY` environment variable → Use legacy if set
2. Detect repository state:
   - `UNINITIALIZED` (no repo) → Use multi-revision (default for new)
   - `LEGACY_ONLY` (only `.git/`) → Use legacy (preserve existing)
   - `BARE_ONLY` (only bare repo) → Use multi-revision
   - `HYBRID` (both exist) → Prefer multi-revision

**Backward compatibility guarantees:**

- Existing repositories continue to work without changes
- `AssetManager` API remains identical
- CLI commands work with both strategies transparently
- Tests pass with minimal modifications
- No forced migration; users opt-in naturally when creating new repos

### Hybrid State Handling

The system gracefully handles hybrid states where both legacy and multi-revision repositories coexist:

- **Detection**: `RepositoryStatus` enum represents all possible states
- **Fallback logic**: Multi-revision strategy can fall back to legacy repo for operations if needed
- **No conflicts**: Strategies are designed to coexist; operations target different directories
- **Explicit control**: Users can force a specific strategy via `setStrategyType()`

### Migration Path

Users naturally migrate as they pull new revisions:

1. **Existing users**: Continue with legacy repos (`LEGACY_ONLY` state detected)
2. **New users**: Get multi-revision by default (`UNINITIALIZED` → multi-revision)
3. **Opt-in migration**: Delete project directory to switch to multi-revision or pull with --migrate
4. **Opt-out**: Set `NXF_SCM_LEGACY=true` to force legacy mode

### Implementation Details

**Key classes:**

- `RepositoryStrategy`: Interface defining repository operations
- `AbstractRepositoryStrategy`: Base class with shared helper methods
- `LegacyRepositoryStrategy`: Direct clone implementation (original behavior)
- `MultiRevisionRepositoryStrategy`: Bare repo + shared clones implementation

**Critical methods:**

- `download()`: Equivalent for both strategies (legacy pulls, multi-revision creates shared clone)
- `getLocalPath()`: Returns appropriate working directory based on strategy
- `getGit()`: Returns appropriate Git instance (legacy git, bare git, or commit git)

### Performance Characteristics

**Disk usage:**
- Legacy: ~100% per repository (full clone with all git objects) + Worktree
- Multi-revision: ~100% for bare + ~100K (.git with alternates) per revision + Worktree per revision

**Operation speed:**
- First download: Similar (both clone from remote)
- Additional revisions: Multi-revision faster (only fetches new objects once, creates cheap clones)
- Switching revisions: Multi-revision instant (different directories), legacy requires checkout

### Known Limitations

- No automatic migration of legacy repositories
- Bare repository overhead even for users who only need one revision
- JGit alternates slightly more complex than worktrees
- Manual cleanup required for old revision clones

## Links
- [GitHub Issue #2870 - Multiple revisions of the same pipeline for concurrent execution](https://github.com/nextflow-io/nextflow/issues/2870)
- [PR #6620 - Implementation of multiple revisions without revisions map](https://github.com/nextflow-io/nextflow/pull/6620)
- Related PRs implementing the multi-revision approach (linked in #6620)

