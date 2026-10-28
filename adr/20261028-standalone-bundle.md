# Standalone bundle with a Java runtime and AOT cache

- Authors: Ben Sherman, Phil Ewels
- Status: draft
- Deciders: TBD
- Date: 2026-10-28
- Tags: launcher, install, java, distribution, performance

Technical Story: [#7696](https://github.com/nextflow-io/nextflow/pull/7696), supersedes [#7697](https://github.com/nextflow-io/nextflow/pull/7697)

## Summary

Each release publishes a standalone bundle per platform, which contains the launcher, a `jlink` Java 25 runtime, the jar and a JVM AOT cache trained in CI. The installed `nextflow` script hands off to the bundle of the selected `NXF_VER`, downloading it instead of the jar when it is available. Nextflow then starts faster, and no longer depends on the Java installed by the user.

## Problem Statement

Nextflow runs on whatever Java it finds, which causes two problems:

- **Installation.** Users must install a compatible Java first. A Java upgrade can break a working install: Homebrew had to pin `openjdk@25` after JDK 27 failed with `Unsupported class file major version 71`.
- **Startup.** Most of the startup of a short run is the JVM loading and linking ~6,000 classes. Java 25 can cache this work ahead of time (JEPs 483, 514, 515), but the cache is only valid for the exact runtime, jar and JVM options it was created with, so it can't be shipped for an arbitrary user Java.

`NXF_VER` adds a constraint: the launcher can run any version from 24.07.0-edge. A runtime must match the selected version, not the launcher. Versions 24.10 and 25.04 fail on Java 25 (`Unsupported class file major version 69`, from Groovy's ASM), while 25.10 and later work. A future runtime upgrade will break older versions in the same way.

## Goals or Decision Drivers

- Each version runs with the runtime and JVM options it was built and tested with, including versions selected with `NXF_VER`.
- A launcher installed today can run versions released later, even if they need a newer Java.
- The AOT cache works on first use, including CI jobs and containers with a fresh `NXF_HOME`.
- Nothing gets worse when a bundle is not available (older versions, other platforms, offline): the launcher falls back to the jar and the system Java.

## Non-goals

- A native binary (GraalVM native image): Nextflow compiles scripts and config files at runtime, and pipelines ship arbitrary Groovy code in `lib/`.
- Windows, Intel macOS and musl Linux bundles (they keep the jar and the system Java).
- Where the bundles are hosted (see open questions).

## Considered Options

1. System Java with an AOT cache created on the user's machine (the original #7696).
2. A Java runtime downloaded by the launcher for `NXF_VER` (#7697).
3. A single runtime bundled with the launcher, used for every `NXF_VER`.
4. A bundle per version and platform, to which the installed launcher hands off.

## Pros and Cons of the Options

### System Java with a local AOT cache

The first `run` records a training profile and the next launch creates the cache in the background.

- Good, because it needs no new release artifacts.
- Bad, because it only applies to users who have Java 25 or later.
- Bad, because it only records in an interactive terminal, so CI jobs and containers with a fresh `NXF_HOME` never get a cache.
- Bad, because the cache is silently rejected on hosts with 128 GB of RAM or more when it was created on a smaller host with a shared `NXF_HOME` (see below).

### Launcher-managed runtime

The launcher downloads a pinned Temurin JRE for `NXF_VER`, from a table of Java versions in the launcher.

- Good, because each version gets a Java it supports, and runtimes are shared between versions.
- Bad, because an installed launcher only knows the versions released before it: it can't pick the runtime of a later version that needs a newer Java.
- Bad, because the JVM options still come from the installed launcher rather than from the version that runs.
- Bad, because the AOT cache must still be created on the user's machine.

### A single bundled runtime

- Good, because it is the simplest layout.
- Bad, because it can't run older versions (24.10 and 25.04 fail on Java 25), nor a later version that needs a newer Java.

### Per-version bundle with a hand-off

- Good, because the runtime, the JVM options and the AOT cache are always those the version was built with, whatever launcher is installed.
- Good, because a launcher installed today can run a later version that needs a newer Java: it downloads that version's bundle, which brings its own runtime.
- Good, because the AOT cache can be created in CI, since the runtime and jar are fixed, so it works on first use.
- Good, because the bundle is self-contained: it can be extracted and run directly, e.g. in container images or air-gapped sites.
- Bad, because there are two layers: the installed script must keep supporting the bundle URL and layout for as long as old launchers exist.
- Bad, because the installed script keeps the jar and system Java path for versions and platforms without a bundle.
- Bad, because each version is a larger download (~85-97 MB instead of 49 MB for the jar), and runtimes are not shared between versions.

## Solution or decision outcome

Option 4: a bundle per version and platform (linux-x64, linux-arm64, macos-arm64), to which the installed launcher hands off.

## Rationale & discussion

### Bundle

```
nextflow-<version>-<os>-<arch>/
  bin/nextflow                  launcher (NXF_VER defaults to the bundled version)
  runtime/                      jlink image of Corretto 25
  lib/nextflow-<version>-one.jar
  lib/nextflow.aot              AOT cache
```

- **Runtime.** It contains the modules found by `jdeps` in the jar, plus those that `jdeps` can't see: service providers (`jdk.zipfs` for Groovy, `jdk.charsets`, `jdk.localedata`, `jdk.naming.dns`), modules used by the core plugins (`jdk.net` for the AWS SDK HTTP client, `jdk.httpserver` for the Seqera Platform and Azure login callbacks), and `jdk.jdwp.agent` for `-remote-debug`. Without `jdk.net`, nf-amazon can't read from S3 at all. The runtime is 124 MB (vs 392 MB for the full JDK). Corretto is used because Temurin 25 ships no `jmods`.
- **AOT cache.** It is trained by running a small pipeline (`bundle/training`) with the bundle launcher, so that it is created with the same JVM options as when it is used. It can only be created by the target runtime, so each platform is built on a native CI runner.
- **Size.** 80-97 MB per platform (tar.gz), of which the runtime is 42 MB, the jar 47 MB and the AOT cache ~12 MB.

### jlink or a stock runtime

The startup gains come from shipping a fixed runtime with a prebuilt AOT cache, which would work the same with a stock JRE: `jlink` itself does not change startup (within noise of a full JDK), nor the per-version runtime. What it provides is size, and control over what the runtime contains:

| runtime (linux-x64) | tar.gz |
|---|---|
| `jlink` image (this bundle) | 42 MB |
| Temurin 25 JRE | 58 MB (53-57 MB on other platforms) |
| Temurin 25 JDK | 134 MB |

The cost is that the module list must be maintained, and a missing module only fails at runtime, when a class of it is used. The basic list (`jdeps` + `jdk.zipfs`) broke S3 in nf-amazon (`NoClassDefFoundError: jdk/net/Sockets`). Third-party plugins run in the same JVM, so they can need a module that no core code or plugin uses. The vendor does not matter: either Corretto or Temurin works, but `jlink` needs the `jmods`, which Temurin 25 does not ship.

`jlink` is kept for the smaller download. A stock JRE remains the fallback if the module list becomes a burden, and the rest of the design does not depend on it.

### Trimming the runtime (follow-up)

`jlink` makes the cost of each module visible, which is an opportunity to review the dependency tree. Compressed size of each module of the runtime, measured by linking it without that module:

| module | tar.gz | used for |
|---|---|---|
| `java.base` | 18 MB | required |
| `jdk.localedata` | 9.2 MB | formatting with non-English locales |
| `java.desktop` native libs | 2.5 MB | never loaded: Groovy only uses `java.beans` (pure Java); excluded by the language server bundle |
| `jdk.charsets` | 1.0 MB | non-standard charsets |
| `jdk.management` | 0.4 MB | `jdeps` (management beans) |
| `java.net.http` | 0.4 MB | `jdeps` (HTTP client) |
| `java.security.jgss` | 0.3 MB | `jdeps` (Kerberos) |
| `jdk.jdwp.agent` | 0.2 MB | `-remote-debug` |
| others | < 0.1 MB each | |

- `jdk.localedata` is most of the extra size. `jlink --include-locales` can keep a subset, if the locales that Nextflow formats with (dates, numbers) are known.
- The `java.desktop` native libraries can be excluded as in the language server.
- `-remote-debug` costs little in size (`jdk.jdwp.agent`), but the feature itself could be reviewed.
- Removing dependencies from the jar (e.g. via `jdeps` per dependency) would shrink both the jar (47 MB) and the module list.

### Launcher

- When the script is in a bundle that contains the jar for `NXF_VER`, it uses the bundled runtime, jar and AOT cache.
- Otherwise, on a supported platform, the installed script hands off to the bundle of `NXF_VER` in `$NXF_HOME/framework/<version>/`. It downloads the bundle when neither the bundle nor the jar of that version is installed, i.e. when it would download the jar today, and falls back to the jar when the bundle is not available. A version without a bundle costs one failed request, once.
- The bundle is extracted into a temporary directory and installed with an atomic symlink, so that concurrent first launches install it once (tested with 4 concurrent launches). Its checksum is verified.
- The hand-off happens after `-download` and `self-update` are handled by the installed script.
- `NXF_JAVA_HOME` still selects another Java, `NXF_BUNDLE=false` disables the hand-off, and `NXF_JVM_ARGS=-XX:AOTMode=off` disables the AOT cache.

### AOT cache validity

The JVM ignores the cache silently (the launcher sets `-Xlog:aot=off`) when:

- the jar mtime differs from training. The bundle archive keeps the file timestamps, which Gradle resets by default for reproducible archives.
- the module options differ, e.g. an extra `--add-opens` in `NXF_OPTS`.
- compressed oops are disabled, which the JVM does by default for a heap of 32 GB or more, i.e. on hosts with 128 GB of RAM or more. The launcher forces `-XX:+UseCompressedOops` unless `NXF_OPTS` sets the heap size.
- ZGC (until JDK 26) or compact object headers are used.

G1 or Parallel GC instead of Serial, the C2 compiler, `-Xmx` up to 32 GB, `-D` properties and moving the bundle to another directory are all fine.

### Results

Measured on linux-x64 (16 CPUs), mean of 5 runs, the same jar and launcher options (C1 compiler, serial GC):

| | full JDK | bundle runtime | bundle runtime + AOT cache |
|---|---|---|---|
| `nextflow info` | 648 ms | 718 ms | 459 ms |
| hello | 1,725 ms | 1,788 ms | 945 ms |
| pipeline with 10 tasks | 1,546 ms | 1,642 ms | 987 ms |
| nf-core/rnaseq `--help` | 8,582 ms | 8,375 ms | 7,553 ms |

The cache saves a roughly constant 0.5-0.8 s per launch, and does not change memory usage or the throughput of long runs. On large pipelines, startup is dominated by script compilation (~5.4 s for rnaseq), which the cache does not cover.

### Distribution and packaging

- **The `-dist` release.** It is a single file with the launcher and the jar, which needs a system Java. It is still needed for platforms without a bundle, and Homebrew builds it from source with its own `openjdk`. It should be kept for now. The bundle could replace it for offline installs on supported platforms, once the bundles are hosted next to the other release files.
- **bioconda** installs the launcher and runs `nextflow -download` at build time, with `NXF_DIST` set inside the package. With the hand-off, this would download the bundle of the build platform into the package, which would break the package on other platforms. The recipe must set `NXF_BUNDLE=false` before the first release with bundles, and the same applies to any package that installs the launcher and pre-downloads the jar.
- **Docker image.** It could use the bundle of its platform instead of installing a JDK.

### Open questions

- Hosting: the workflow attaches the bundles to the GitHub release, and the launcher downloads them from there (`NXF_BUNDLE_BASE` overrides the location). They could be uploaded next to the jar under `www.nextflow.io/releases` instead.
- `self-update` run from an extracted bundle replaces the bundle's launcher with the plain one. It keeps working, through the hand-off, but the bundle is no longer self-contained.
- Old versions accumulate in `$NXF_HOME/framework` as today; `self-update -clean` removes them, including their bundles.
- The training pipeline covers a typical run. Commands such as `lint` load classes that it does not cover, which are loaded as usual.

## Links

- Supersedes [#7697](https://github.com/nextflow-io/nextflow/pull/7697) (launcher-managed runtime) and the local AOT cache of the original [#7696](https://github.com/nextflow-io/nextflow/pull/7696).
- [nextflow-io/language-server#191](https://github.com/nextflow-io/language-server/pull/191): the same bundle for the language server, whose AOT cache is created on the user's machine.
- [JEP 483](https://openjdk.org/jeps/483), [JEP 514](https://openjdk.org/jeps/514), [JEP 515](https://openjdk.org/jeps/515).
