# Plugin Lockfile Integrity

- Authors: Claude <noreply@anthropic.com>
- Status: draft
- Date: 2026-07-27
- Tags: plugins, security, registry, integrity

## Summary

Introduce a committed lockfile (`plugins.lock`) that pins the SHA-512 hash of each resolved plugin **archive**, computed over the actual archive bytes. The lockfile is populated automatically as a side effect of the normal plugin download — the first time a coordinate is fetched it is pinned (trust-on-first-use) — exactly like `go.sum` / `package-lock.json`, with **no dedicated command**. The trusted hash lives in the pipeline repository under version control and is not re-fetched from the registry on every run. During plugin resolution Nextflow re-hashes the downloaded — or retained — archive and verifies it against the locked hash.

Verification is network-free **whenever the archive being verified is present locally** (always true on a cold download; true on a warm cache only if the archive was retained). This ADR is explicit about the one case where an archive is not present — a cache extracted before this feature existed — and specifies an offline-safe behavior for it that never triggers a download (see *Migration and first adoption* and *Verification flow*).

## Problem Statement

Nextflow resolves plugins from a registry (or mirror), downloads a zip archive, and extracts it into the local plugin cache under `$id-$version/`. Today there is no persistent, user-owned record of what a plugin archive is *supposed* to hash to. The registry returns a `sha512sum` during resolution, but that value is fetched fresh on every run and is only ever compared against the download in transit — so it can detect corruption on the wire, never a compromise of the source of truth itself.

**What is already covered today.** In-transit archive integrity is *not* an open gap: `HttpPluginRepository` wires pf4j's `CompoundVerifier` (`HttpPluginRepository.groovy:134-135`), which checks the downloaded archive against the registry-supplied `sha512sum` (`:222`) at download time. A corrupted or MITM-mangled download that no longer matches the registry's own hash is already rejected. The lockfile does **not** claim novelty there.

The genuine, non-redundant gaps the lockfile closes are:

- **Runtime trust in the registry's own hash**: today the *expected* value is whatever the registry asserts on this run. A compromised registry (or mirror) that serves a tampered archive together with a matching tampered hash passes the existing `CompoundVerifier` check — because both sides of the comparison come from the same untrusted source. The lockfile removes the registry from the trust path at run time.
- **Silent version drift**: a coordinate that resolves to a different artifact than it did last week produces no signal today.
- **Reproducibility / committed baseline for teams that vendor or mirror plugins**: there is no local, authoritative, VCS-reviewed record of what each coordinate must hash to, independent of the registry.

The core issue is that the *trusted* hash must be owned by the pipeline and committed to version control, so verification depends on the repository — which the team controls and reviews — rather than on the registry being honest at run time.

## Goals or Decision Drivers

- **Local, network-free verification when the archive is present**: once the lockfile exists, checking a locally-present archive against the locked hash requires no network access. Network is used only to *fetch* a plugin that is not yet cached, exactly as today. This ADR does **not** claim verification is unconditionally offline — see the migration case below.
- **Content-addressed trust anchor**: the pinned hash is computed over the downloaded archive bytes at first download, so the baseline is a genuine content hash (the `go.sum` / `package-lock.json` model), not a re-recording of a registry claim.
- **User-owned trust anchor**: the authoritative hash is committed to the pipeline repository and reviewed through normal VCS workflows (pull requests, code review, blame).
- **Reproducibility**: a given commit of a pipeline resolves to exactly the same plugin archives, or fails loudly.
- **Zero breakage by default**: pipelines without a lockfile behave exactly as they do today. The feature is opt-in by the mere presence of the file, and its default enforcement level mirrors PR #7308's warn-first rollout philosophy (see *Modes and rollout*).
- **Never turn verification into a network dependency**: the verification path must never re-fetch from the registry as a remedy — that would re-introduce the exact runtime dependency this design rejects in Option A.
- **Honest scope**: the mechanism must defend what it can actually defend (the archive, and only as a re-extraction gate on a warm cache) and not overclaim protection it does not provide (the extracted tree that actually runs).

## Non-goals

- **Extracted-tree hashing**: hashing or re-verifying the contents of the extracted `$id-$version/` directory is out of scope. That surface is covered by the `PluginSecurity` directory guard in PR #7308 (see *Threat model* and *Known residual* below).
- **Transitive / dependency lock resolution**: only the coordinates declared in the pipeline config are locked. There is no dependency-graph resolution or locking of plugins pulled in indirectly.
- **Registry-server changes**: this ADR is entirely client-side. The registry contract is unchanged.
- **Signature or provenance verification**: the lockfile pins a hash, not a cryptographic signature or attestation. Provenance metadata (the `url`) is recorded for auditability only.
- **Catching in-transit download corruption**: already handled by pf4j's `CompoundVerifier` against the registry-supplied `sha512sum`. Not a job the lockfile duplicates.

## Considered Options

### Option A: Load-time checksum comparison against the registry

Fetch the plugin, then re-query the registry for the expected `sha512sum` at load time and compare. No committed lockfile.

- Good, because there is nothing new to commit or maintain in the pipeline repo.
- Bad, because it makes the **registry a runtime dependency** of every plugin load — breaking air-gapped and offline execution, the exact environments where integrity matters most.
- Bad, because the registry becomes both the source of the artifact *and* the source of the trusted hash: a compromised registry can serve a tampered archive with a matching tampered hash and defeat the check entirely (this is exactly the `CompoundVerifier` behavior that already exists today).
- Bad, because it provides no reproducibility anchor — the "expected" value can change out from under a pipeline between runs with no committed record.

### Option B: Extracted directory tree-hash

Compute a hash over the extracted `$id-$version/` directory tree and pin that.

- Good, because it would detect edits to already-extracted files (cache poisoning of the unpacked plugin) — the surface the archive hash cannot see.
- Bad, because tree hashing is heavier: it must walk and hash the full extracted tree, and is sensitive to extraction non-determinism (file ordering, timestamps, permissions, symlinks) across platforms and unzip implementations.
- Bad, because it duplicates the protection the `PluginSecurity` directory guard in PR #7308 already provides for the extracted tree.
- Bad, because it would tempt an implementation to re-extract and re-hash on every run, adding cost to the warm-cache hot path.
- Deferred: the extracted-tree surface is real and is what actually executes, but it is owned by the #7308 guard, not by the lockfile. (Note: a variant of Option B — deriving the lock hash from the extracted tree — is the only way to make verification offline on a pre-existing cache with no retained archive; that is called out where relevant below but remains deferred.)

### Option C: Committed archive-hash lockfile with download-gate and retain-and-re-verify (adopted)

Commit a `plugins.lock` file pinning the SHA-512 of each resolved plugin **archive**, hashed from the archive bytes at generation time. On a cold cache, gate the download against the locked hash. Retain the verified archive in the cache and re-verify it on warm-cache runs.

- Good, because the trusted hash is owned by the pipeline and committed to VCS — the registry is not trusted at verification time.
- Good, because the pinned hash is content-addressed (computed over the bytes), so it survives registry compromise on all runs *after* the lock was honestly generated (trust-on-first-use).
- Good, because it adds a reproducibility / drift anchor the registry cannot silently move.
- Good, because verification of a locally-present archive is network-free.
- Good, because archive hashing is deterministic and cheap (one hash of one zip), with no extraction-ordering pitfalls.
- Bad (accepted), because it requires **retaining** the archive in the cache — a new persisted artifact with size and lifecycle cost (see *Tradeoff: retaining the archive*), and one an attacker can simply ignore.
- Bad (accepted), because on a warm cache the thing that actually executes is the extracted tree, not the retained zip — so re-hashing the zip protects only a *future* re-extraction, not the current run (see *Threat model* and *Known residual*). Warm-run integrity of the executed code rests on the #7308 guard.
- Bad (accepted), because a cache extracted *before* this feature has no retained archive, so its archive cannot be verified offline; the design degrades safely rather than re-downloading (see *Migration and first adoption*).

## Solution or decision outcome

**Option C — committed archive-hash lockfile with download-gate and retain-and-re-verify** — is the recommended approach. It places a content-addressed trust anchor in the pipeline repository, keeps verification network-free whenever the archive is locally present, never re-fetches from the registry as a verification remedy, and composes cleanly with the extracted-directory guard from PR #7308 rather than duplicating it.

## Rationale & discussion

### Threat model

The lockfile's genuine contribution (beyond what pf4j's `CompoundVerifier` already does at download time) is:

- **Removing runtime trust in the registry-supplied hash** — a compromised registry/mirror that serves a tampered archive with a matching tampered hash passes today's `CompoundVerifier`, but fails against the committed content hash. This holds on every run after the lock was honestly generated.
- **Silent version drift** — a coordinate resolving to a different artifact fails against the committed hash.
- **A committed reproducibility baseline** for teams that vendor or mirror plugins.

What the lockfile does **not** defend:

- **In-transit corruption** — already caught by `CompoundVerifier`; not a lockfile novelty.
- **A poisoned *extracted* cache** — this is the important honesty point. On a warm cache Nextflow loads and executes the extracted `$id-$version/` tree; the retained zip is not what runs. An attacker who can write to the plugin cache will simply modify the extracted tree and leave the retained zip untouched, and the lockfile check passes while poisoned code executes. Re-hashing the archive therefore provides **essentially zero protection for what actually runs on a warm cache**; it only detects tampering of the archive itself, which matters solely if that archive is later re-extracted. Warm-run integrity of the executed plugin rests **entirely** on the `PluginSecurity` directory guard in **PR #7308**, which is a required complement, not optional.

The two features compose:

- the **lockfile** guarantees the *archive* you obtained (and retained) is the archive you committed to, and gates the cold-cache download and any future re-extraction;
- the **#7308 directory guard** guarantees the *extracted tree* — the code that actually executes — has not been altered after extraction.

This split is stated plainly: the lockfile must not be marketed as "surviving cache poisoning." It survives archive tampering; the directory guard survives extracted-file tampering.

### Lockfile format and location

- **File name**: `plugins.lock`, in the pipeline project root next to `nextflow.config`, committed to VCS.
- **Format**: JSON, chosen for diff-friendliness under code review.
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

- `sha512` is the hash of the plugin **archive** (the zip), **computed over the downloaded archive bytes**. It is a content hash the pipeline owns, not a re-recording of the registry's `sha512sum`. (No `url` field is stored: it is not available at the point the archive is retained, and it would be provenance-only — never a trust input.)

### Generation — automatic on first download (trust-on-first-use)

There is **no dedicated `lock` command**. The lockfile is populated as a side effect of the normal plugin download, exactly as `go` writes `go.sum` and `npm` writes `package-lock.json` on first install:

- The feature is dormant unless a `plugins.lock` file is present. To start, create an empty one (`touch plugins.lock`) and run the pipeline once.
- On a **cold-cache download** the archive is being fetched anyway; Nextflow hashes the retained bytes. If the coordinate is **missing** from the lock, its computed hash is **appended** (trust-on-first-use). Reviewing the resulting diff and committing it establishes the baseline.
- Because the pinned value is computed from the bytes actually received (not copied from registry metadata), it is a genuine content hash — the `go.sum` model — which is what makes the "survives registry compromise" property honest.
- Trust model: TOFU. Pinning trusts the registry at the moment a coordinate is first seen; every run thereafter is anchored to the committed content hash and no longer trusts the registry. This is why the pinned lockfile is meant to be reviewed and committed to VCS.

Auto-pinning only ever **adds a missing coordinate**. An entry already present is **never silently rewritten** — a downloaded archive whose hash differs from a committed entry is a verification *failure* (mode-gated below), not a silent update, exactly as `go.sum` refuses to quietly change a recorded hash. Re-pinning a legitimately changed plugin is an explicit action: delete the stale entry and re-run.

### Verification flow (in `PluginUpdater`)

Verification happens during plugin resolution in `PluginUpdater`. **No branch of this flow ever re-fetches from the registry as a remedy** — fetching happens only to obtain a plugin that is genuinely absent from the cache, exactly as today.

- **Cold cache (download path)**: the archive is downloaded to obtain the plugin regardless of the lockfile. After fetching the zip, compute its SHA-512. If the coordinate is **already locked**, compare against the committed entry — the **lock**, not the registry-supplied `sha512sum`, is authoritative. If the coordinate is **not yet locked**, append it (trust-on-first-use). Retain the zip in the cache (`$id-$version.zip` next to the extracted `$id-$version/` directory) so later runs can re-verify it offline.
- **Warm cache, archive retained (extracted dir + retained zip present)**: re-hash the retained zip against the lock, offline. This gates a *future* re-extraction only; see the threat model for why it does not protect the currently-executing extracted tree.
- **Warm cache, archive absent (extracted dir present, no retained zip)** — the universal state for every cache extracted before this feature (see *Migration and first adoption*): the plugin is already present and functional, so **do not download anything**. Verification of the archive is simply not possible offline for this cache. Behavior is:
  - `strict` / `warn`: emit a one-time notice that the archive is unavailable for lockfile verification and that a fresh download (or re-vendoring the archive) is needed to establish an offline-verifiable baseline. **Do not abort** and **do not re-download** — a present, functioning plugin is never failed solely because its archive is missing and cannot be fetched offline. Integrity of the extracted tree for this run is provided by the #7308 directory guard.
  - `off`: skip silently.
- **Coordinate not in the lock**: on a cold download it is **auto-pinned** (trust-on-first-use, see *Generation*); on a warm cache with no retained archive there are no bytes to pin, so it is a silent no-op.

Verification of a **locally-present** archive requires no network. The one case where the archive is not present (a pre-feature cache) is handled without any network access, by design — it never falls back to a download.

#### Known residual (documented, not overclaimed)

Re-hashing the retained zip proves the **archive** is intact. On a warm run it does **not** protect the code that actually executes: Nextflow runs the already-extracted `$id-$version/` tree, and re-hashing the zip verifies an artifact that is not the thing being run. For a warm run the archive re-verify is therefore effectively inert for the executed code — its only value is gating a subsequent re-extraction. Detection of extracted-file tampering — the integrity of what actually runs — is the responsibility of the `PluginSecurity` directory guard in PR #7308. This residual is expected and by design.

#### Tradeoff: retaining the archive

Today Nextflow unzips the archive and immediately deletes it (`PluginUpdater.groovy:267-269`); a warm cache returns the extracted dir with no archive (`:259-262`). This ADR changes that for locked plugins by keeping `$id-$version.zip` alongside the extracted directory. Costs and caveats:

- **Cache size**: roughly doubles on-disk footprint per plugin (compressed archive + extracted tree). Acceptable for the reproducibility/offline benefit; could be scoped to locked plugins only.
- **Lifecycle**: the retained zip must be cleaned up with the plugin directory and re-written on re-download.
- **Attacker can ignore it**: retaining the zip adds no protection for the running code — an attacker edits the extracted tree and leaves the zip untouched, and the check still passes. This is precisely why the #7308 directory guard is a required complement.

### Migration and first adoption

Because current code deletes the archive right after extraction (`PluginUpdater.groovy:267-269`) and warm-cache resolution returns the extracted directory with no archive (`:259-262`), **every plugin cache that predates this feature has no retained archive**. This is the *universal initial state* on first adoption, not an edge case.

Consequences, stated plainly:

- On the first lockfile-enabled run against a pre-existing cache — including an air-gapped one — plugins are already extracted and functional. The archive-absent branch above applies: Nextflow does **not** re-download, does **not** abort in strict mode, and simply notes that an offline-verifiable archive baseline has not yet been established. No network is required and no false abort occurs.
- An offline-verifiable archive baseline is established the next time each plugin is downloaded through the gate (cold cache) once a `plugins.lock` file exists, at which point the archive is retained and the coordinate pinned.
- **Air-gapped teams that vendor or mirror plugins** typically ship the extracted `$id-$version/` directories, not the `.zip` archives. To get *archive-level* offline verification in such an environment, the archives must be vendored too — i.e. the mirror/vendor step must include the retained `$id-$version.zip` files (produced by running the pipeline once against the registry in a connected environment, then committing/shipping the archives alongside the lockfile). Absent that, air-gapped runs fall into the archive-absent branch and rely on the #7308 directory guard for the extracted tree — which is safe and network-free, but is not archive verification. The only alternative that would make archive-free caches offline-verifiable is deriving the lock hash from the extracted tree (deferred Option B).

### Modes and rollout

- **Opt-in by presence**: if no `plugins.lock` exists, the feature is dormant and behavior is unchanged — zero breakage for existing pipelines.
- When the file exists, behavior is controlled by `NXF_PLUGINS_LOCK_MODE`, reusing the `warn`/`strict`/`off` tri-state plumbing introduced by PR #7308's `NXF_PLUGINS_STRICT_MODE`, but **independent** from it (the two knobs can be set separately).

Mode gates only the **hash-mismatch** outcome. A coordinate *missing* from the lock is auto-pinned (see *Generation*), not gated; an *absent* archive is never gated (never aborts).

| Mode | Hash mismatch (locked entry) | Coordinate missing from lock | Archive absent (pre-feature cache) |
|------|------------------------------|------------------------------|-------------------------------------|
| `strict` | **Abort** with a re-pin hint | Auto-pin (cold) / no-op (warm) | Notice only, proceed (never abort) |
| `warn` (default when the file is present) | Log a warning once and proceed | Auto-pin (cold) / no-op (warm) | Notice only, proceed |
| `off` | Skip verification | Auto-pin (cold) / no-op (warm) | Skip silently |

**Default is `warn` when a lockfile is present**, deliberately mirroring PR #7308's warn-first rollout philosophy (`NXF_PLUGINS_STRICT_MODE` defaults to `warn`). This was a considered choice: an earlier draft proposed `strict`-when-present on the reasoning that a committed lockfile expresses intent to enforce. That was rejected because it diverges from the #7308 rollout philosophy this feature otherwise mirrors, and because it produces a surprising hard failure in the one gated case — a plugin whose committed lock entry is legitimately stale (e.g. the archive was re-released) would hard-abort until the entry is deleted and re-pinned. `warn`-by-default surfaces the drift without breaking the run; teams that want enforcement opt into `strict` explicitly (a staged `warn` → `strict` rollout), exactly as with #7308. (Note the auto-pin model already removes the most common friction: a *new or bumped* coordinate is pinned on first download rather than aborting, even in `strict`.)

Independently of mode, a **present and functioning plugin is never aborted solely because its archive is missing** and cannot be fetched offline (the archive-absent column above) — this avoids a false-positive failure on precisely the air-gapped caches the feature is meant to serve.

### Reused and new components

| Component | Module | Change |
|-----------|--------|--------|
| `PluginLockFile` (new) | nf-commons | Read/write/round-trip of `plugins.lock`; blank/malformed-file handling |
| `PluginLockVerifier` (new) | nf-commons | Auto-pin on first download (TOFU); re-verify retained archive; archive-absent handling; mode gating |
| `PluginUpdater` | nf-commons | Cold-cache: retain zip + verify/pin; warm-cache: re-verify retained zip. No dedicated command |
| `HttpPluginRepository` | nf-commons | No change — `CompoundVerifier` still checks downloads against the registry `sha512sum` at fetch time |
| `NXF_PLUGINS_LOCK_MODE` | nf-commons | New env var, tri-state (`strict`/`warn`/`off`), default `warn`, independent of `NXF_PLUGINS_STRICT_MODE` |

### Relationship to PR #7308

PR #7308 introduces the `PluginSecurity` directory guard, which protects the **extracted** plugin directory — the code that actually executes — and the `NXF_PLUGINS_STRICT_MODE` tri-state plumbing (default `warn`). This ADR's lockfile is the **sibling** feature protecting the **archive**, reusing the same mode-plumbing pattern and the same warn-first default under a separate, independent switch. Together they cover both integrity surfaces — archive and extracted tree — without either overclaiming the other's protection. Critically, warm-run integrity of executing code depends on the #7308 guard; the lockfile does not substitute for it.

## Testing

- **`PluginLockFile`**: read, write, and round-trip of `plugins.lock`; blank/`touch`ed file parses as empty; malformed file throws.
- **Auto-pin (TOFU)**: an enabled but empty lock, given a downloaded archive for an unlocked coordinate, appends the **byte-computed** hash to the file on disk; a dormant (no file) verifier never pins.
- **Verification**:
  - matching hash passes;
  - mismatch **aborts** in `strict`, **warns once** in `warn`, is **ignored** in `off`;
  - **archive-absent (pre-feature cache)**: with a locked coordinate but no retained zip, the run proceeds in `strict` with a notice, performs **no download**, and never aborts.
- **No-network guarantee**: verification and pinning operate purely on locally-present bytes; no branch re-fetches from the registry as a remedy.
