# Plugin Lockfile Integrity

- Authors: Paolo Di Tommaso <paolo.ditommaso@gmail.com>
- Status: proposed
- Date: 2026-07-27
- Tags: plugins, security, registry, integrity

## Summary

Add a committed lockfile, `plugins.lock`, that pins a canonical SHA-512 hash of each plugin's extracted directory, the unpacked `$id-$version/` tree that Nextflow loads and executes. Nextflow populates the lockfile the first time it sees a plugin, the way `go` writes `go.sum` and `npm` writes `package-lock.json`. There is no command to run. The pinned hash lives in the pipeline repository under version control, and Nextflow does not re-fetch it from the registry on every run. Before loading a plugin, Nextflow re-hashes its extracted directory and compares the result against the pinned hash.

Hashing the code that runs, rather than the archive it came from, covers two problems with one check: a tampered download or a compromised registry, and a poisoned cache directory, where a lower-trust user edits already-extracted files on a shared filesystem. Verification is offline. It reads no ownership or permission bits, so it stays quiet on a legitimate shared cache.

## Problem Statement

Nextflow resolves plugins from a registry or mirror, downloads a zip archive, extracts it into the local plugin cache under `$id-$version/`, and loads the extracted code. That leaves two gaps.

1. **Cache poisoning of the executed code.** When a `$id-$version/` directory already exists, `PluginUpdater.load0()` skips the download and loads the extracted tree as it stands. On a shared, writable plugin cache, where `NXF_PLUGINS_DIR` points at a group-writable or world-writable location, a lower-trust local user can pre-populate or edit an official pinned coordinate such as `nf-amazon-3.10.0/` with their own class and jar files. Those files then execute inside the victim's Nextflow process.
2. **Runtime trust in the registry.** The registry returns a `sha512sum` during resolution, but Nextflow fetches that value fresh on every run and only compares it against the download in transit, through pf4j's `CompoundVerifier` (`HttpPluginRepository.groovy:134-135` and `:222`). A compromised registry or mirror that serves a tampered archive together with the tampered archive's hash passes the check, and silent version drift produces no signal at all. Nothing on the client records what a plugin is *supposed* to be.

What is already covered: `CompoundVerifier` catches in-transit corruption against the registry-supplied hash, and the lockfile does not duplicate that. The lockfile covers what happens after the download, namely the extracted code, and it removes runtime trust in the registry's own hash.

The baseline has to be owned by the pipeline and committed to version control, and it has to describe the artifact that executes, the extracted tree. Verification then depends on the repository, which the team controls and reviews, instead of on the registry being honest at run time or the cache being untouched.

## Goals or Decision Drivers

- **Cover the executed code.** The check must protect the extracted `$id-$version/` tree that Nextflow loads, not only the archive it was unpacked from.
- **One check, both problems.** The same content hash should answer registry integrity and cache poisoning, rather than shipping two features.
- **No false positives.** Verification must not read directory ownership or permissions. A legitimate shared or admin-managed cache must produce no output at all.
- **Local verification.** Once the lockfile exists, checking a plugin needs no network access. The registry is never required at verification time.
- **A pinned hash the user owns.** The hash is computed over the files themselves and committed to the pipeline repository, reviewed like any other change.
- **No breakage by default.** Pipelines without a lockfile behave as they do today. The presence of the file turns the feature on, and a mismatch warns by default, matching PR #7308's warn-first rollout.

## Non-goals

- **A dedicated `lock` command.** Nextflow populates the lockfile on first load. There is no `nextflow plugins lock` subcommand to run or maintain.
- **Transitive dependency locking.** Only coordinates that are loaded get pinned. Nothing resolves a dependency graph.
- **Registry-server changes.** The work is entirely client-side and the registry contract is unchanged.
- **Signature or provenance verification.** The lockfile pins a content hash, not a signature or an attestation.
- **Catching in-transit download corruption.** pf4j's `CompoundVerifier` already handles it.

## Considered Options

### Option A: Load-time checksum comparison against the registry

Fetch the plugin, re-query the registry for the expected `sha512sum` at load time, and compare. No committed lockfile.

- Bad, because every plugin load then needs the registry, which breaks air-gapped and offline execution. Those are the environments where integrity matters most.
- Bad, because the registry supplies both the artifact and the hash it is judged by. A compromised registry defeats the check.
- Bad, because it records nothing for reproducibility and does nothing about cache poisoning.

### Option B: Committed archive-hash lockfile

Commit a `plugins.lock` pinning the SHA-512 of each plugin archive, the zip, verified on download and, by retaining the zip, on warm-cache runs.

- Good, because it is a committed, offline, content-addressed record for registry integrity and drift.
- Bad, because on a warm cache the archive is not what executes. Nextflow runs the extracted tree, so re-hashing a retained zip cannot detect edits to the unpacked files. It does nothing about cache poisoning and would still need PR #7308's directory guard alongside it.
- Bad, because it means keeping the archive next to the extracted tree, roughly doubling the cache footprint. It also introduces an "archive absent" case, covering every cache that predates the feature and every admin cache that ships only extracted directories, which produces notices without adding protection.

### Option C: Committed extracted-tree-hash lockfile (adopted)

Commit a `plugins.lock` pinning a canonical SHA-512 over the contents of the extracted `$id-$version/` directory. Before loading a plugin, re-hash its directory and compare against the committed entry.

- Good, because it hashes the bytes that execute, so one hash covers registry drift and cache poisoning together.
- Good, because it reads no ownership or permission bits. It is silent on any legitimate cache, whether private, shared read-only, or admin service-account, and speaks only on a real content mismatch. That removes the warning noise a directory-ownership guard produces.
- Good, because it keeps no archive, so the cache does not double in size, and there is no "archive absent" case: a plugin that loads always has an extracted tree.
- Good, because verifying a directory that is already on disk needs no network.
- Bad, accepted: it re-hashes the extracted tree. For large cloud plugins that cost is real, which motivates the caching optimisation described below.
- Bad, accepted: like any lockfile, the first pin trusts the artifact present when the coordinate is first seen.

## Solution

Option C, the committed extracted-tree-hash lockfile, is adopted. It pins the code that runs, keeps verification offline, produces no ownership-based false positives, and covers both Option B's goal and PR #7308's goal with one check. It supersedes the archive-hash approach and the directory-ownership guard.

## Rationale & discussion

### Threat model

- **Cache poisoning of the executed code.** A lower-trust user edits the extracted `$id-$version/` tree on a shared, writable cache. Re-hashing the tree at load detects any change to the files, whoever owns the directory and whatever its permissions. It therefore catches what an ownership heuristic misses, such as a group-writable cache owned by a trusted service account, and it never mis-fires on a legitimate read-only shared cache.
- **Registry or mirror compromise, and drift.** A tampered archive extracts to a different tree, so its hash no longer matches the committed entry. That holds on every run after the lock was generated honestly.
- **In-transit corruption.** `CompoundVerifier` already catches it. The lockfile adds nothing here.

The first pin of each coordinate is trust-on-first-use. It trusts the tree present when the coordinate is first seen, which on a cold download was just fetched and checked in transit against the registry hash. Once the pin is committed, every later run compares against that reviewed baseline and trusts neither the registry nor the cache.

### Lockfile format and location

- **File name.** `plugins.lock`, in the pipeline project root next to `nextflow.config`, committed to version control. The project directory is the one holding the main script, which for a pipeline pulled from a Git repository is its local clone under `$NXF_HOME/assets`. `CmdRun` passes it to the plugin system once the script has been resolved. It is not the launch directory, so a lockfile committed by one pipeline can never apply to another.
- **Concurrency.** Nextflow writes the file to a sibling temp file and moves it into place, so a concurrent reader never sees a truncated file. Two runs pinning the same missing coordinate at once is last-writer-wins, and the next run repairs it, since both compute the same hash.
- **Format.** JSON, for readable diffs under code review, with a sorted key order.
- **Keying.** By resolved `id@version`.

```json
{
  "version": 1,
  "plugins": {
    "nf-amazon@2.0.0": {
      "sha512": "cbc4..."
    }
  }
}
```

- `sha512` is a canonical hash of the extracted directory tree. For every regular file, in sorted relative-path order, the digest absorbs the relative path and the file bytes, but not timestamps or permissions. The hash is therefore stable across extractions and platforms while still catching any change to the code.

### Generation: automatic on first load

There is no command to run. Nextflow populates the lockfile as a side effect of loading a plugin, the way `go` writes `go.sum` and `npm` writes `package-lock.json`.

- The feature is dormant unless a `plugins.lock` file is present. To start, create an empty one with `touch plugins.lock` and run the pipeline once.
- The first time a coordinate is loaded and it is missing from the lock, Nextflow appends its extracted-tree hash. Review the resulting diff and commit it to establish the baseline.
- The pinned value comes from the files on disk, so it is a real content hash.

Auto-pinning only adds a coordinate that is missing. Nextflow never rewrites an entry that is already there. A plugin whose tree differs from a committed entry is a verification failure, gated by the mode described below, not a silent update, in the same way `go.sum` refuses to quietly change a recorded hash. Re-pinning a plugin that changed legitimately is deliberate: delete the stale entry and run again.

### Verification flow (in `PluginUpdater`)

`load0()` verifies each plugin once, immediately before `loadPluginFromPath()`, so a fresh download and a reused cache take the same path. No branch re-fetches from the registry as a remedy. Nextflow reaches the network only to obtain a plugin that is genuinely absent from the cache, exactly as it does today.

- Compute the canonical hash of the extracted `$id-$version/` directory.
- If the coordinate is not in the lock, append it (trust-on-first-use).
- If it matches the committed entry, proceed.
- If it differs, the mode decides (below).

Verifying a directory that is already on disk needs no network. There is no "artifact absent" case, because a plugin that is loading has an extracted directory by definition.

### Performance and the re-hash cost

When a lockfile is present, every load re-hashes the extracted directory. For small plugins that is nothing next to JVM and Nextflow start-up. For large cloud plugins, which bundle SDKs of tens to hundreds of MB, it costs real time.

The first implementation hashes in full every time and does not pre-optimise. A follow-up, added only once the cost is shown to matter, skips the re-hash when nothing changed:

- Keep a local, uncommitted sidecar file, for example under `$NXF_PLUGINS_DIR`, recording per plugin the last verified tree hash and a cheap fingerprint of the tree, the set of `(relative-path, size, mtime)`.
- On load, stat the files without reading them and rebuild the fingerprint. If both the fingerprint and the lock entry are unchanged, reuse the last verified result and skip the full hash.
- Anyone who can write the files can forge `mtime` and `size`, so trust the fingerprint only when the cache cannot change underneath you, that is, a directory you own that is not writable by group or others. On a shared or writable cache, ignore the fingerprint and always hash in full. Ownership therefore picks the strategy, fast or full, and never produces a warning, so the noise problem of an ownership guard does not come back.

Steady-state launches on a private cache stay nearly free, the full hash runs wherever tampering is possible, and nothing is printed except on a real mismatch.

### Modes and rollout

- No `plugins.lock`, no behaviour change: the feature is dormant.
- When the file exists, `NXF_PLUGINS_LOCK_MODE` decides what a mismatch does. It reuses the `warn`/`strict`/`off` tri-state of PR #7308's `NXF_PLUGINS_STRICT_MODE` but is independent of it.

| Mode | Tree hash mismatch (locked entry) | Coordinate missing from lock |
|------|-----------------------------------|------------------------------|
| `strict` | **Abort** with a re-pin hint | **Abort**, an unpinned coordinate is refused, not pinned |
| `warn` (default when the file is present) | Log a warning once and proceed | Auto-pin |
| `off` | Feature disabled: no hashing, no pinning, no reporting | |

In `warn`, Nextflow auto-pins a coordinate that is missing from the lock, and the mode gates only the mismatch. `strict` gates both. It refuses a coordinate that reaches the loader with no lock entry, because an unpinned coordinate can resolve past the committed version to a tree that is already extracted in a shared cache, where the download short-circuit means pf4j's verifier never ran either, and auto-pinning it would silently record attacker-supplied bytes. Populate the lock in `warn` mode, commit it, then switch to `strict`. `off` is a complete way out, since Nextflow does not even read the lock file, so a run can always be unblocked without deleting a file committed in the pipeline repository. The default is `warn` for a gentle rollout, and `strict` is there for teams that want the run to fail closed. Unlike a directory-ownership guard, both defaults stay silent in normal operation: a mismatch is a rare, actionable event, not per-run noise on a healthy cache.

### Relationship to PR #7308

PR #7308 proposed a `PluginSecurity` directory guard that flags a plugin directory as untrusted when it is foreign-owned or world-writable, to defend the extracted tree on shared caches. The extracted-tree hash defends the same thing by content instead of by ownership. That is stronger, because it catches any modification whoever owns the directory, including the group-writable, trusted-owner case the heuristic misses, and quieter, because it never warns on a legitimate shared or read-only cache. It also covers registry integrity, which the guard did not. This ADR therefore supersedes PR #7308. Keeping the plugin cache private, one per user, remains the zero-config advice, and the lockfile is the committed, enforceable layer for teams that want it.

### Reused and new components

| Component | Module | Change |
|-----------|--------|--------|
| `PluginLockFile` (new) | nf-commons | Read/write/round-trip of `plugins.lock`; blank and malformed file handling |
| `PluginLockVerifier` (new) | nf-commons | Canonical extracted-tree hash; auto-pin on first load; mismatch mode gating |
| `PluginUpdater` | nf-commons | Verify or pin the extracted directory in `load0()` before loading. No retained archive, no dedicated command |
| `HttpPluginRepository` | nf-commons | No change. `CompoundVerifier` still checks downloads against the registry `sha512sum` at fetch time |
| `NXF_PLUGINS_LOCK_MODE` | nf-commons | New env var, tri-state (`strict`/`warn`/`off`), default `warn`, independent of `NXF_PLUGINS_STRICT_MODE` |

## Testing

- **`PluginLockFile`**: read, write, round-trip; a blank or `touch`ed file parses as empty; a malformed file and a future format version both abort with an `AbortOperationException`; writes are atomic.
- **`sha512Tree`**: identical trees in different locations hash equal; any change to a file changes the hash; the output is a 128-char lowercase hex digest.
- **Auto-pin**: given an enabled but empty lock and an extracted directory for an unlocked coordinate, the tree hash is appended to the file on disk; a dormant verifier, with no file present, never pins.
- **Verification**:
  - a matching tree passes;
  - a mismatch aborts in `strict` and warns once in `warn`, and `off` leaves the verifier dormant;
  - a coordinate missing from the lock aborts in `strict` and is auto-pinned in `warn`;
  - in `warn`, pinning a coordinate whose plugin id is already locked at a different version logs a version-drift warning;
  - a lock file that cannot be written, for example in a read-only checkout, warns instead of aborting the run;
  - cache poisoning: pinning a good tree and then editing an extracted file in place is detected as a mismatch, whoever owns the directory.
- **No network**: verification and pinning read only local files, and no branch re-fetches from the registry.
