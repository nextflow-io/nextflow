# Nextflow-managed Java runtime

- Authors: Phil Ewels
- Status: draft
- Deciders: TBD
- Date: 2026-09-26
- Tags: launcher, install, java, distribution

Technical Story: [#2028](https://github.com/nextflow-io/nextflow/issues/2028), [#4380](https://github.com/nextflow-io/nextflow/issues/4380)

## Summary

Nextflow needs users to install a compatible Java before they can run it. That is the most common source of install problems. With this change, the `nextflow` launcher downloads a pinned, known-good Java runtime on first use and uses it by default. When it can't, it falls back to the system Java. Package-managed installs (conda, Homebrew, container images, HPC modules) keep providing their own Java.

## Problem Statement

Today the launcher runs whatever Java it finds, from `NXF_JAVA_HOME`, `JAVA_HOME` or `PATH`. Each Nextflow release only works with a range of Java versions, and users have to match it by hand:

- **No Java installed.** New users hit an error before their first run. The install docs spend most of their length on installing Java through SDKMAN.
- **Java too new.** A system or package-manager Java upgrade can break Nextflow overnight. For example, Homebrew had to pin `openjdk@25` after JDK 27 failed with `Unsupported class file major version 71`.
- **Java too old or slow.** Old distro JDKs miss startup and GC improvements. Old `NXF_VER` releases need an older Java than the one the user has installed.
- **Unreliable Java packaging.** Some channels ship a problematic Java (see the conda discussion in [#4380](https://github.com/nextflow-io/nextflow/issues/4380)).

## Goals or Decision Drivers

- **It must always work.** A fresh machine with no Java runs Nextflow. When the managed runtime can't be used, behaviour is never worse than today.
- A known-good Java for each Nextflow version, including the old versions selected with `NXF_VER`.
- `self-update` and `NXF_VER` keep working.
- Package managers keep full control of versions and dependencies.
- Negligible per-launch overhead once installed.
- Mirrors and air-gapped installs are supported.

## Non-goals

- Removing the JVM dependency (for example GraalVM native-image). Runtime compilation of scripts, `lib/` Groovy code and JVM plugins rule this out.
- Changing the downstream recipes (bioconda, Homebrew). This ADR only defines the contract they rely on.
- Garbage collection of old runtimes in `$NXF_HOME/jre`.

## Considered Options

1. System Java only (status quo), with better docs.
2. Bundled runtime: a per-platform tarball with a `jlink`-trimmed runtime next to the `-dist` launcher.
3. Launcher-managed runtime: the launcher downloads a pinned JRE on first use.
4. Native binary (GraalVM native-image).

## Pros and Cons of the Options

### System Java only

- Good, because nothing changes.
- Bad, because every failure mode above remains, and the install docs stay SDKMAN-heavy.

### Bundled runtime

A prototype built a 17-module `jlink` runtime. The bundle was about 77 MB per platform versus 49 MB for today's `-dist`. Startup was within 5% of a full JDK. HTTPS and plugin downloads worked.

- Good, because nothing is downloaded at runtime, which suits air-gapped sites.
- Bad, because the jar is embedded in the script, so `NXF_VER` is silently ignored.
- Bad, because `self-update` replaces the script with the plain launcher, which then can't find Java.
- Bad, because a trimmed module list can miss modules that third-party plugins need.
- Bad, because every release needs about 4 extra assets of about 77 MB each.

### Launcher-managed runtime

- Good, because it removes the Java prerequisite for direct installs without changing the distribution format.
- Good, because each Nextflow version gets a matching Java, including old `NXF_VER` releases.
- Good, because `self-update` and `NXF_VER` keep working: the runtime lives in `$NXF_HOME`, the same way framework jars do.
- Good, because it uses the full Temurin JRE, so no plugin can fail on a missing module.
- Bad, because the first launch downloads about 40 to 60 MB and depends on the download host being reachable.
- Bad, because it changes behaviour: users who set `JAVA_HOME` no longer get that Java by default.

### Native binary

- Good, because it would give instant startup and a single file.
- Bad, because Nextflow compiles pipeline scripts and config files to bytecode at runtime, and pipelines ship arbitrary Groovy in `lib/`. A closed-world native image can't run either, so this would mean replacing the Groovy runtime.

## Solution or decision outcome

The launcher manages the Java runtime by default for direct installs (option 3). Package-managed installs opt out and keep supplying their own Java. The bundled tarball (option 2) stays available as an optional channel for air-gapped sites.

## Rationale & discussion

### Java resolution order

1. `NXF_JAVA_HOME`, if set.
2. The managed runtime, unless `NXF_JAVA_MANAGED=false`.
3. The existing lookup (`JAVA_HOME`, `/usr/libexec/java_home`, `PATH`), only as a fallback. The launcher prints a warning that says why the managed runtime was not used and names both overrides.

With a managed runtime, the launcher skips its own Java-version check, because the runtime is pinned to a known-good version.

### Compatibility matrix

The Java major version is chosen from `NXF_VER`. The table is derived from the Java versions each release's launcher accepted (`version_check`), picking the newest LTS in range:

| `NXF_VER` | Managed Java |
|---|---|
| ≥ 25.09.0-edge | 25 LTS |
| 23.09.2-edge – 25.04.x | 21 LTS |
| 21.10.0 – 23.09.1-edge | 17 LTS |
| 18.10.0 – 21.04.x | 11 LTS |
| < 18.10 | none (system Java) |

The table lives in the launcher. It pins exact Temurin releases (for example `25.0.4.1+1`), so a given launcher always installs the same runtime. An `NXF_VER` newer than the launcher knows about gets the newest row. Each release updates the pins as part of `make releaseInfo`. "Known good" should be enforced in CI by running each supported release line against its pinned runtime.

### Download and verification

- **Source:** Eclipse Temurin JRE tarballs, from the stable GitHub release URLs (`github.com/adoptium/temurin<N>-binaries`). This covers glibc and musl Linux on x64 and aarch64, and macOS on x64 and aarch64. Other platforms fall back to the system Java.
- **Verification:** the tarball's `.sha256.txt` must match before extraction. The extracted `bin/java -version` must run before the install is published. Today the Nextflow jar itself is downloaded without any checksum.
- **Mirrors:** `NXF_JAVA_BASE` overrides the host, with the same path layout; `file://` also works. Before general availability, the runtimes should be mirrored under the Nextflow release host. Firewalls that allow the Nextflow jar download would then also allow the runtime.

### Install layout and concurrency

- **Layout:** each runtime is extracted into a temporary directory under `$NXF_HOME/jre/`. The launcher then publishes it by atomically creating the symlink `<version>-<os>-<arch>[-musl]`.
- **Concurrency:** the symlink avoids the `mv dir existing_dir` nesting trap. When concurrent first launches race, one wins. The others delete their temporary copy and use the winner's.
- **Failure cleanup:** interrupted, corrupt or non-running downloads are never published, and an EXIT trap removes their temporary directories.
- **Permissions:** directories are world-readable, so a shared `NXF_HOME` works for every user.

### TLS trust

- **The problem:** the Temurin runtime ships its own `cacerts`. On its own, it would reject corporate proxy CAs that are installed in the OS trust store.
- **Linux:** the launcher passes `-Djavax.net.ssl.trustStore` pointing at the distro's Java trust store (`/etc/pki/java/cacerts` or `/etc/ssl/certs/java/cacerts`) when that file exists. This matches what the distro's own OpenJDK does. A trust store the user sets in `NXF_OPTS` still wins.
- **macOS:** `-Djavax.net.ssl.trustStoreType=KeychainStore-ROOT` would pick up CAs installed by device management (MDM). That store type only exists in JDK 23+, so it needs a version guard and testing on managed Macs first. It is not in the prototype.
- **Not covered:** users who imported CAs into their own JDK's `cacerts` need `NXF_JAVA_MANAGED=false`, or a trust store set in `NXF_OPTS`.

### Package-managed installs

The package manager owns the Nextflow version, the Java dependency and updates:

| Channel | Java | Nextflow version | `self-update` / `NXF_VER` |
|---|---|---|---|
| Direct install (`get.nextflow.io`) | managed | `NXF_VER`, `self-update` | work |
| bioconda, Homebrew, Docker image, HPC modules, Spack/Nix | package's own (`NXF_JAVA_MANAGED=false`) | fixed by the package | should be disabled with a pointer to the package manager |
| Air-gapped tarball (optional) | bundled (`jlink`) | fixed | disabled |

Proposed, not yet prototyped: packagers set a single variable, for example `NXF_INSTALLER=conda`. It implies `NXF_JAVA_MANAGED=false`, and it makes `self-update` print the right upgrade command instead of overwriting a file the package manager owns.

### Cost

Measured with the runtime installed and the version cache warm:

- **Per-launch overhead:** 89.4 ms vs 89.5 ms for the launcher alone, and 381 ms vs 383 ms for `-version`. This is within noise.
- **First launch:** about 3 s extra on a fast connection.

Runtime tarball sizes per platform:

| JRE | Size |
|---|---|
| 25 | 42–62 MB |
| 21 | 42–52 MB |
| 17 | 38–47 MB |
| 11 | 38–44 MB |

### Prototype test results

All tests passed in bare containers (Debian bookworm-slim, Alpine, Rocky 9) and on macOS:

| Scenario | Result |
|---|---|
| No Java: `-version`, `run nextflow-io/hello`, second launch | managed runtime installed once, second launch silent |
| Alpine/musl (x64 and aarch64), busybox `wget` only | works |
| `JAVA_HOME` set to another JDK | managed runtime used; `NXF_JAVA_MANAGED=false` restores the system JDK |
| No network, system Java present / absent | warning plus fallback / clear error naming the overrides |
| Read-only `NXF_HOME` | warning plus fallback |
| Bad checksum, truncated tarball, missing checksum file, runtime that does not run | not installed, falls back, `jre/` stays empty |
| 8 concurrent first launches on one `NXF_HOME` | all succeed, exactly one install |
| `NXF_VER` from 18.10.1 to 25.10.4 (one or more per matrix row) | correct runtime selected, `-version` and `info` succeed |
| `self-update`, then launch | same runtime, no new download |
| Private CA installed with `update-ca-certificates` | HTTPS works through the OS trust store |

### Known limitations and open questions

- **Concurrent first launches:** HPC job arrays or CI fan-out download the runtime N times. The results are correct, but the downloads could hit rate limits. Running one launch on a login node first avoids it. A lock would add stale-lock risk.
- **Offline nodes:** they retry the download on every launch before falling back. That takes up to 10 s when a firewall drops packets silently. A negative cache could avoid it.
- **Old runtimes:** they accumulate in `$NXF_HOME/jre` when the pins move.
- **Stdout pollution:** the existing jar download message (`get()`) goes to stdout, so the first run of `nextflow config` still shows it. The runtime download message goes to stderr.
- **`JAVA_CMD` in tasks:** the exported `JAVA_CMD` points at the managed runtime, which task environments inherit.
- **Hosting:** mirroring the runtimes under the Nextflow release host is recommended but not yet decided.
- **macOS keychain trust:** see TLS trust above.

## Links

- Related: [#2028](https://github.com/nextflow-io/nextflow/issues/2028) native package distribution, [#2951](https://github.com/nextflow-io/nextflow/issues/2951) standalone core distribution (done), [#4380](https://github.com/nextflow-io/nextflow/issues/4380) installation improvements.
- Related launcher changes, proposed separately: fix the Java-version cache never being written, and a JVM AOT cache for faster startup on Java 25+. A managed Java 25 runtime makes the AOT cache apply by default.
