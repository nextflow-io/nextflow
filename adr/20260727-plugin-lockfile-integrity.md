# Plugin Lockfile Integrity

- Authors: Claude <noreply@anthropic.com>
- Status: draft
- Date: 2026-07-27
- Tags: plugins, security, registry, integrity

## Summary

Introduce a committed lockfile (`plugins.lock`) that pins a canonical SHA-512 hash of each plugin's **extracted directory** — the unpacked `$id-$version/` tree that Nextflow actually loads and executes. The lockfile is populated automatically the first time a plugin is seen (trust-on-first-use), exactly like `go.sum` / `package-lock.json`, with **no dedicated command**. The trusted hash lives in the pipeline repository under version control and is not re-fetched from the registry on every run. Before loading a plugin, Nextflow re-hashes its extracted directory and verifies it against the committed hash.

Because it hashes the code that runs rather than the archive it came from, a single mechanism covers **both** integrity surfaces: a tampered/compromised download or registry, **and** a poisoned cache directory (a lower-trust user editing already-extracted files on a shared filesystem). Verification is fully offline, and it uses no ownership or permission heuristics, so it does not produce false-positive warnings on legitimate shared caches.

## Problem Statement

Nextflow resolves plugins from a registry (or mirror), downloads a zip archive, extracts it into the local plugin cache under `$id-$version/`, and loads the extracted code. Two integrity gaps exist:

1. **Cache poisoning (the executed code).** When a `$id-$version/` directory already exists, `PluginUpdater.load0()` short-circuits the download entirely and loads the extracted tree as-is. On a shared, writable plugin cache (`NXF_PLUGINS_DIR` pointed at a group-writable or world-writable location), a lower-trust local user can pre-populate or edit an official pinned coordinate (e.g. `nf-amazon-3.10.0/`) with attacker-controlled class/jar files that then execute inside the victim's Nextflow process.
2. **Runtime trust in the registry.** The registry returns a `sha512sum` during resolution, but that value is fetched fresh every run and only ever compared against the download in transit (pf4j's `CompoundVerifier`, `HttpPluginRepository.groovy:134-135` / `:222`). A compromised registry/mirror that serves a tampered archive with a matching tampered hash passes that check, and silent version drift produces no signal. There is no persistent, user-owned record of what a plugin is *supposed* to be.

**What is already covered today.** In-transit corruption of a download is caught by `CompoundVerifier` against the registry-supplied hash. The lockfile does not duplicate that. Its job is the *post-download* surface: the extracted code, plus removing runtime trust in the registry's own hash.

The core issue is that the *trusted* baseline must be owned by the pipeline and committed to version control, and it must describe the artifact that actually executes — the extracted tree — so verification depends on the repository (which the team controls and reviews) rather than on the registry being honest at run time or on the cache being un-tampered.

## Goals or Decision Drivers

- **Cover the executed code.** The check must protect the extracted `$id-$version/` tree that Nextflow loads, not merely the archive it was unpacked from.
- **One mechanism, both surfaces.** Supply-chain/registry integrity and cache-poisoning integrity should be answered by the same content hash, not two features.
- **No false positives, no noise.** Verification must not rely on directory ownership or permission heuristics; a legitimate shared or admin-managed cache must be silent. Output happens only on a genuine mismatch.
- **Local, network-free verification.** Once the lockfile exists, checking a plugin requires no network access; the registry is never a runtime dependency of verification.
- **Content-addressed, user-owned trust anchor.** The pinned hash is computed over the actual files and committed to the pipeline repository, reviewed through normal VCS workflows.
- **Zero breakage by default.** Pipelines without a lockfile behave exactly as today; the feature is opt-in by the mere presence of the file, defaulting to `warn` on mismatch (mirroring PR #7308's warn-first rollout).

## Non-goals

- **A dedicated `lock` command.** The lockfile is populated automatically on first load; there is no `nextflow plugins lock` subcommand to run or maintain.
- **Transitive / dependency lock resolution.** Only coordinates that are actually loaded are pinned; there is no dependency-graph resolution.
- **Registry-server changes.** Entirely client-side; the registry contract is unchanged.
- **Signature / provenance verification.** The lockfile pins a content hash, not a cryptographic signature or attestation.
- **Catching in-transit download corruption.** Already handled by pf4j's `CompoundVerifier`.

## Considered Options

### Option A: Load-time checksum comparison against the registry

Fetch the plugin, then re-query the registry for the expected `sha512sum` at load time and compare. No committed lockfile.

- Bad, because it makes the **registry a runtime dependency** of every plugin load — breaking air-gapped and offline execution, the environments where integrity matters most.
- Bad, because the registry becomes both the source of the artifact *and* the source of the trusted hash: a compromised registry defeats the check entirely.
- Bad, because it provides no reproducibility anchor and does nothing for cache poisoning.

### Option B: Committed archive-hash lockfile

Commit a `plugins.lock` pinning the SHA-512 of each plugin **archive** (the zip), verified on download and (by retaining the zip) on warm-cache runs.

- Good, because it is a committed, offline, content-addressed anchor for supply-chain integrity and drift.
- Bad, because on a warm cache **the archive is not what executes** — Nextflow runs the extracted tree. Re-hashing a retained zip cannot detect edits to already-unpacked files, so it does **not** cover cache poisoning (the executed code) at all. It would still need PR #7308's directory guard as a complement.
- Bad, because it requires **retaining the archive** next to the extracted tree — roughly doubling cache footprint — and introduces an "archive absent" case (every pre-feature cache, and admin caches that ship only extracted dirs) that produces notices without adding protection.

### Option C: Committed extracted-tree-hash lockfile (adopted)

Commit a `plugins.lock` pinning a canonical SHA-512 over the **extracted `$id-$version/` directory contents**. Before loading a plugin, re-hash its directory and compare to the committed entry.

- Good, because it hashes exactly the bytes that get executed, so a **single** hash covers supply-chain/drift **and** cache poisoning.
- Good, because it uses **no ownership/permission heuristic** — it is silent on any legitimate cache (private, shared read-only, admin service-account) and speaks only on a real content mismatch. This eliminates the false-positive/warning-noise problem of a directory-ownership guard.
- Good, because it needs **no retained archive** (no doubled footprint) and has **no "archive absent" case** — the extracted tree is always present when a plugin loads.
- Good, because verification of a locally-present directory is network-free.
- Bad (accepted), because it must re-hash the extracted tree; for large cloud plugins this cost is non-trivial and motivates the caching optimisation described below.
- Bad (accepted), because, like any lockfile, the first pin is trust-on-first-use — it trusts the artifact present when the coordinate is first seen.

## Solution or decision outcome

**Option C — committed extracted-tree-hash lockfile** — is adopted. It places a content-addressed, user-owned trust anchor over the code that actually runs, keeps verification network-free, produces no ownership-based false positives, and with one mechanism covers both the supply-chain surface (Option B's goal) and the cache-poisoning surface (PR #7308's goal). It therefore **supersedes** both the archive-hash approach and the directory-ownership guard.

## Rationale & discussion

### Threat model

- **Cache poisoning of the executed code.** A lower-trust user edits the extracted `$id-$version/` tree on a shared/writable cache. Re-hashing the tree at load detects any change to the files, **independently of who owns the directory or its permissions** — so it catches cases an ownership heuristic misses (e.g. a group-writable cache owned by a trusted service account) and never mis-fires on a legitimate read-only shared cache.
- **Registry/mirror compromise and drift.** A tampered archive extracts to a different tree, so its hash no longer matches the committed entry. This holds on every run after the lock was honestly generated.
- **In-transit corruption.** Already caught by `CompoundVerifier`; not a lockfile novelty.

The first pin of each coordinate is **trust-on-first-use**: it trusts the tree present when first seen (which, on a cold download, was just fetched and checked in transit against the registry hash). Once committed to VCS, every later run is anchored to that reviewed baseline and no longer trusts the registry or the cache.

### Lockfile format and location

- **File name**: `plugins.lock`, in the pipeline project root next to `nextflow.config`, committed to VCS.
- **Format**: JSON, chosen for diff-friendliness under code review, with a stable (sorted) key order.
- **Keying**: by resolved `id@version`.

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

- `sha512` is a canonical hash of the **extracted directory tree**: for every regular file, in sorted relative-path order, the digest absorbs the relative path and the file bytes (not timestamps or permissions), so it is stable across extractions and platforms while detecting any change to executable content.

### Generation — automatic on first load (trust-on-first-use)

There is **no dedicated command**. The lockfile is populated as a side effect of loading a plugin, as `go` writes `go.sum` and `npm` writes `package-lock.json`:

- The feature is dormant unless a `plugins.lock` file is present. To start, create an empty one (`touch plugins.lock`) and run the pipeline once.
- The first time a coordinate is loaded and it is **missing** from the lock, its extracted-tree hash is **appended**. Reviewing the resulting diff and committing it establishes the baseline.
- The pinned value is computed from the files on disk, so it is a genuine content hash.

Auto-pinning only ever **adds a missing coordinate**. An entry already present is **never silently rewritten** — a plugin whose tree differs from a committed entry is a verification *failure* (mode-gated below), not a silent update, exactly as `go.sum` refuses to quietly change a recorded hash. Re-pinning a legitimately changed plugin is explicit: delete the stale entry and re-run.

### Verification flow (in `PluginUpdater`)

Verification happens once per plugin in `load0()`, immediately before `loadPluginFromPath()`, covering both a fresh download and a reused (warm) cache with the same code path. **No branch ever re-fetches from the registry as a remedy** — network access happens only to obtain a plugin genuinely absent from the cache, exactly as today.

- Compute the canonical hash of the extracted `$id-$version/` directory.
- If the coordinate is **not in the lock** → append it (trust-on-first-use).
- If it **matches** the committed entry → proceed.
- If it **differs** → mode-gated (below).

Verification of a locally-present directory requires no network. There is no "artifact absent" case: if a plugin is being loaded, its extracted directory exists by definition.

### Performance and the re-hash cost

Verification re-hashes the extracted directory on every load when a lockfile is present. For small plugins this is negligible next to JVM and Nextflow start-up; for large cloud plugins (bundled SDKs of tens to hundreds of MB) it is material.

The initial implementation performs the full hash each time and does **not** pre-optimise. A follow-up optimisation, added only once the cost is shown to matter, avoids the re-hash when nothing changed:

- Keep a **local, non-committed** sidecar (e.g. under `$NXF_PLUGINS_DIR`) recording, per plugin, the last verified tree hash and a cheap fingerprint of the tree — the set of `(relative-path, size, mtime)`.
- On load, stat the files (no reads) and rebuild the fingerprint; if it is unchanged and the lock entry is unchanged, reuse the last verified result and skip the full hash.
- **Safety of the shortcut:** `mtime`/`size` are forgeable by anyone who can write the files, so the fingerprint shortcut is trusted **only when the cache cannot change under you** — i.e. a directory you own that is not writable by group or others. On a shared/writable cache the fingerprint is not trusted and the full hash always runs. Ownership is thus used to select the *strategy* (fast vs full), never to emit a warning — so the noise problem of an ownership guard does not reappear.

This keeps steady-state launches on a private cache nearly free, runs the full hash exactly where tampering is possible, and stays silent in all cases except a genuine mismatch.

### Modes and rollout

- **Opt-in by presence**: no `plugins.lock` → dormant, behaviour unchanged, zero breakage.
- When the file exists, mismatch behaviour is controlled by `NXF_PLUGINS_LOCK_MODE`, reusing the `warn`/`strict`/`off` tri-state plumbing pattern from PR #7308's `NXF_PLUGINS_STRICT_MODE`, but independent from it.

| Mode | Tree hash mismatch (locked entry) | Coordinate missing from lock |
|------|-----------------------------------|------------------------------|
| `strict` | **Abort** with a re-pin hint | Auto-pin (trust-on-first-use) |
| `warn` (default when the file is present) | Log a warning once and proceed | Auto-pin |
| `off` | Skip verification | Auto-pin |

Mode gates only the **mismatch** outcome. A coordinate missing from the lock is auto-pinned, not gated. Default is `warn` for a friendly rollout; `strict` is opt-in for teams that want a fail-closed guarantee. Unlike a directory-ownership guard, **either default is silent in normal operation** — a mismatch is a real, rare, actionable event, not per-run noise on healthy caches.

### Relationship to PR #7308

PR #7308 proposed a `PluginSecurity` directory guard that flags a plugin directory as untrusted when it is foreign-owned or world-writable, to defend the extracted tree on shared caches. This ADR's extracted-tree hash defends the **same** surface by content instead of by ownership, which is both stronger (it catches any modification regardless of ownership, including the group-writable-trusted-owner case the heuristic misses) and quieter (it never warns on a legitimate shared/read-only cache). It also covers the supply-chain surface the guard did not. This mechanism therefore **supersedes PR #7308**, which can be closed in its favour. The zero-config baseline for the extracted-tree surface remains the operational guidance to keep the plugin cache private (per-user); the lockfile is the committed, enforceable layer on top for teams that want it.

### Reused and new components

| Component | Module | Change |
|-----------|--------|--------|
| `PluginLockFile` (new) | nf-commons | Read/write/round-trip of `plugins.lock`; blank/malformed-file handling |
| `PluginLockVerifier` (new) | nf-commons | Canonical extracted-tree hash; auto-pin on first load (TOFU); mismatch mode gating |
| `PluginUpdater` | nf-commons | Verify/pin the extracted directory in `load0()` before loading. No retained archive, no dedicated command |
| `HttpPluginRepository` | nf-commons | No change — `CompoundVerifier` still checks downloads against the registry `sha512sum` at fetch time |
| `NXF_PLUGINS_LOCK_MODE` | nf-commons | New env var, tri-state (`strict`/`warn`/`off`), default `warn`, independent of `NXF_PLUGINS_STRICT_MODE` |

## Testing

- **`PluginLockFile`**: read, write, round-trip; blank/`touch`ed file parses as empty; malformed file throws.
- **`sha512Tree`**: identical trees in different locations hash equal; any file-content change changes the hash; output is a 128-char lowercase hex digest.
- **Auto-pin (TOFU)**: an enabled but empty lock, given an extracted directory for an unlocked coordinate, appends the tree hash to the file on disk; a dormant (no file) verifier never pins.
- **Verification**:
  - matching tree passes;
  - mismatch **aborts** in `strict`, **warns once** in `warn`, is **ignored** in `off`;
  - **cache poisoning**: pinning a good tree then editing an extracted file in place is detected as a mismatch — regardless of directory ownership.
- **No-network guarantee**: verification and pinning operate purely on local files; no branch re-fetches from the registry.
