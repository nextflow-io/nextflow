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
- Bad, because the AOT cache is only valid for the exact JVM build it was created with (e.g. 25.0.1, not 25.0.2), which conflicts with sharing runtimes: either the runtime build is pinned per version, which gives a runtime per version stored separately from the jar, or the cache is created on the user's machine.
- Bad, because a cache created on the user's machine makes the first run slower, and each install runs an artifact that was never tested.
- Bad, because the launcher must manage runtime versions and remove unused ones, and the runtime and jar are separate downloads that can fail independently.

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
- Bad, because each version is a larger download (80-101 MB instead of 49 MB for the jar), and runtimes are not shared between versions.

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
- **Size.** 80-101 MB per platform (tar.gz). On linux-x64, the runtime is 43 MB, the jar 44 MB and the AOT cache 14 MB.

### jlink or a stock runtime

The startup gains come from shipping a fixed runtime with a prebuilt AOT cache, which would work the same with a stock JRE: `jlink` itself does not change startup (within noise of a full JDK), nor the per-version runtime. What it provides is size, and control over what the runtime contains:

| runtime (linux-x64) | tar.gz |
|---|---|
| `jlink` image (this bundle) | 42 MB |
| Temurin 25 JRE | 58 MB (53-57 MB on other platforms) |
| Temurin 25 JDK | 134 MB |

The cost is that the module list must be maintained, and a missing module only fails at runtime, when a class of it is used. The basic list (`jdeps` + `jdk.zipfs`) broke S3 in nf-amazon (`NoClassDefFoundError: jdk/net/Sockets`). Third-party plugins run in the same JVM, so they can need a module that no core code or plugin uses. The vendor does not matter: either Corretto or Temurin works, but `jlink` needs the `jmods`, which Temurin 25 does not ship.

A stock JRE would be unpacked as `runtime/` in CI, in the same per-version bundle, with the same hand-off and AOT cache, so nothing else changes. Compared with `jlink`:

- Good, because there is no module list to maintain, and no module can be missing at runtime.
- Good, because it needs no `jmods`, so any vendor's JRE can be used.
- Neutral, because it ships the default CDS archive of the JDK, which the AOT cache supersedes.
- Bad, because the download is ~16 MB larger per version and platform, for modules Nextflow does not use.

`jlink` is kept for the smaller download. A stock JRE remains the fallback if the module list becomes a burden, and the rest of the design does not depend on it.

### Trimming the bundle (follow-up)

`jlink` makes the cost of each module visible, which is an opportunity to review the dependency tree. Measured on the linux-x64 bundle of 26.10.3 (compressed sizes):

| part | current | floor |
|---|---|---|
| runtime | 43.3 MB | 29.1 MB with the modules that core and Groovy need; 18.7 MB for `java.base` alone |
| jar | 43.6 MB | ~14 MB for the Nextflow modules and Groovy |
| AOT cache | 14.2 MB | smaller, in proportion to the classes loaded during training |
| total | 101 MB | ~55 MB |

**Runtime.** Size saved by linking the runtime without each module:

| module | saved | used for |
|---|---|---|
| `jdk.localedata` | 9.6 MB | formatting with non-English locales; not referenced by the jar |
| `java.desktop` | 9.0 MB | required for `java.beans` (Groovy, snakeyaml); its native libs (2.6 MB) are never loaded |
| `jdk.charsets` | 1.0 MB | non-standard charsets |
| `java.net.http` | 0.4 MB | required (lib-httpx, npr-client) |
| `java.security.jgss` | 0.3 MB | Kerberos in jgit |
| `jdk.jdwp.agent` | 0.2 MB | `-remote-debug` |
| others | < 0.1 MB each | |

Applied in steps: excluding the `java.desktop` native libs gives 40.6 MB, adding `--include-locales=en` gives 31.1 MB, and keeping only the modules needed by core and Groovy gives 29.1 MB. `java.base` is mostly the JVM itself (`libjvm.so`, 28 MB raw), and Corretto only ships the server VM.

**Jar.** The Nextflow modules and Groovy are 14 MB of the 43.6 MB jar. The largest dependencies:

| dependency | size | notes |
|---|---|---|
| Bouncy Castle | 12.4 MB | only used by nf-k8s (embedded since 26.08) to read PEM keys, in `SSLUtils.java` |
| jgit | 3.4 MB | `pull` and `clone`; only user of `java.security.jgss` |
| Guava | 3.3 MB | used across nf-commons and nextflow |
| Jackson + json-schema-validator | 3.1 MB | used by one class; overlaps with Gson and groovy-json |
| ANTLR 4 tool + StringTemplate + ANTLR 3 | 1.8 MB | the grammar generator of nf-lang, leaked onto the runtime classpath by the Gradle `antlr` plugin |
| Woodstox | 1.8 MB | registered as the `XMLInputFactory`, which the JDK provides |
| GPars | 1.7 MB | |
| Apache Commons (compress, lang3, io, codec) | 2.7 MB | |
| others (logback, javax.mail, jsoup, kryo, leveldb, ...) | ~5 MB | |

The cheapest reductions:

- Exclude the `java.desktop` native libs, as in the language server, once nf-console is removed (2.6 MB).
- Keep only English locale data with `--include-locales=en` (9.5 MB), if formatting dates and numbers in English for all users is acceptable.
- Remove the ANTLR tool from the runtime classpath (1.8 MB).
- Read PEM keys with the JDK in nf-k8s instead of Bouncy Castle (12.4 MB).

Together they bring the bundle to ~75 MB. Beyond that, each dependency would need to be replaced or removed, which would also shrink the AOT cache. `-remote-debug` costs little in size, but the feature itself could be reviewed.

### Static compilation of core classes (follow-up)

Most core classes are dynamic Groovy: in the `nextflow` module, 658 of 794 classes (excluding closures) contain dynamic call sites. Groovy links each one at runtime through its own `invokedynamic` bootstrap, which the AOT cache does not archive (it only covers JDK lambdas and string concatenation), so every launch relinks them and generates the classes behind them.

Converting the classes that link the most call sites in a hello run (~35 classes, e.g. `TaskProcessor`, `TaskRun`, `ProcessConfig`, the script params, `DAG`, `HistoryFile`) to `@CompileStatic` gave, with the AOT cache:

| | dynamic | `@CompileStatic` |
|---|---|---|
| hello | 1,127-1,158 ms | 1,024-1,041 ms (~10% faster) |
| `nextflow info` | 454-457 ms | 437-446 ms (within noise) |
| call sites linked (hello) | 2,548 | 1,235 |
| hidden classes generated (hello) | 634 | 297 |
| AOT cache | 53.6 MB | 53.1 MB |
| jar | 49.0 MB | 49.0 MB |

- Size barely changes: core classes are ~10 MB of the 49 MB jar (raw), and ~18% of the classes loaded at startup, so neither static compilation nor Java would shrink the bundle much. The Groovy runtime stays in any case, for scripts and config.
- Startup improves, and more is possible: almost all remaining call sites come from statically compiled code, where Groovy 4 still emits dynamic calls for `as` coercions and implicit casts (to `String`, `boolean`, `Map`). Removing all of them would save an estimated ~200 ms on hello. Java avoids them by construction.
- The conversion is not mechanical. `@CompileStatic` changes behaviour silently: property reads on map-like config classes bypass their getters (e.g. `config.cacheable` on `ProcessConfig`), overload resolution and `in` checks differ, and two Groovy compiler bugs (`clone()` on an interface, varargs wrapping) produced broken code. Unit tests caught most of these, but a wider conversion would need a review of the generated bytecode as well.

### Launcher

- When the script is in a bundle that contains the jar for `NXF_VER`, it uses the bundled runtime, jar and AOT cache.
- Otherwise, on a supported platform, the installed script hands off to the bundle of `NXF_VER` in `$NXF_HOME/framework/<version>/`. It downloads the bundle when neither the bundle nor the jar of that version is installed, i.e. when it would download the jar today, and falls back to the jar when the bundle is not available. A version without a bundle costs one failed request, once.
- The bundle is extracted into a temporary directory and installed with an atomic symlink, so that concurrent first launches install it once (tested with 4 concurrent launches). It is only installed if its checksum matches and its runtime runs (it may not, e.g. on a `noexec` mount). The download times out after 10 s without a connection, and reports to stderr, so that the output of commands such as `nextflow config` is not affected.
- The hand-off happens after `-download` and `self-update` are handled by the installed script.
- `NXF_JAVA_HOME` still selects another Java, `NXF_JAVA_BUNDLE=false` disables the hand-off, and `NXF_JVM_ARGS=-XX:AOTMode=off` disables the AOT cache.
- Mirrors and air-gapped installs can set `NXF_JAVA_BUNDLE_BASE` to another host with the same layout (`<base>/v<version>/<bundle>.tar.gz` and `.sha256`), including a `file://` URL.

Changes in behaviour, on supported platforms:

- `JAVA_HOME` and the `java` on the `PATH` are no longer used, only `NXF_JAVA_HOME`.
- The exported `JAVA_CMD`, which tasks of the local executor inherit, points to the bundled runtime. Nextflow itself does not use it.

### TLS trust

The bundled runtime has its own CA certificates, so it would not trust the CAs installed in the OS, e.g. by a corporate proxy that inspects TLS. Unless `NXF_OPTS` or `NXF_JVM_ARGS` set a trust store, the launcher uses the trust store of the OS, as the Java of the OS would:

- **Linux:** the Java trust store of the distro, `/etc/pki/java/cacerts` (RHEL and derivatives, generated by `update-ca-trust`) or `/etc/ssl/certs/java/cacerts` (Debian and Ubuntu, only if `ca-certificates-java` is installed), when it exists. Otherwise the bundled CA certificates are used.
- **macOS:** the keychain (`-Djavax.net.ssl.trustStoreType=KeychainStore-ROOT`, JDK 23+), which includes the root certificates installed by device management. Not yet tested on a managed Mac.

Users who added CAs to the trust store of their own JDK need `NXF_JAVA_HOME`, or a trust store in `NXF_OPTS`.

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
- **bioconda** installs the launcher and runs `nextflow -download` at build time, with `NXF_DIST` set inside the package. With the hand-off, this would download the bundle of the build platform into the package, which would break the package on other platforms. The recipe must set `NXF_JAVA_BUNDLE=false` before the first release with bundles, and the same applies to any package that installs the launcher and pre-downloads the jar.
- **Docker image.** It is based on `amazoncorretto:25-al2023` (glibc), so it could use the bundle of its platform instead of installing a JDK.
- **Package managers.** A package manager owns the Nextflow version, the Java dependency and updates. Instead of each recipe knowing about `NXF_JAVA_BUNDLE=false`, packages could set a single variable, e.g. `NXF_INSTALLER=conda`, which disables the hand-off and makes `self-update` print the upgrade command of the package manager instead of overwriting a file that it owns. Not implemented yet.

### Testing

In bare containers (linux-x64), with the bundles served from a local mirror:

| Scenario | Result |
|---|---|
| Rocky 9, CentOS 7 (glibc 2.17), Debian slim, Ubuntu 24.04, no Java installed | bundle installed on first launch, pipeline runs with the AOT cache (`-XX:AOTMode=on`) |
| Alpine (musl) | no bundle, same error as today when no Java is installed |
| No network, `file://` mirror | works |
| No network, or dropped packets (blackhole IP) | falls back to the jar immediately, or after the 10 s timeout |
| Read-only `NXF_HOME` with the bundle installed | works |
| `NXF_HOME` on a `noexec` mount | bundle not installed, jar and system Java used |
| Private CA added with `update-ca-trust` (Rocky 9) | HTTPS works through the OS trust store; fails with the bundled CA certificates only |
| 4 concurrent first launches | all succeed, one bundle installed |
| `NXF_VER` without a bundle (25.10.9) | jar and system Java used, no error |

CI builds each bundle on its platform, runs the training pipeline from an extracted copy with `-XX:AOTMode=on`, and reads a file over HTTPS.

### Open questions

- Hosting: the workflow attaches the bundles to the GitHub release, and the launcher downloads them from there (`NXF_JAVA_BUNDLE_BASE` overrides the location). They could be uploaded next to the jar under `www.nextflow.io/releases` instead.
- `self-update` run from an extracted bundle replaces the bundle's launcher with the plain one. It keeps working, through the hand-off, but the bundle is no longer self-contained.
- Old versions accumulate in `$NXF_HOME/framework` as today; `self-update -clean` removes them, including their bundles.
- The training pipeline covers a typical run. Commands such as `lint` load classes that it does not cover, which are loaded as usual.

## Links

- Supersedes [#7697](https://github.com/nextflow-io/nextflow/pull/7697) (launcher-managed runtime) and the local AOT cache of the original [#7696](https://github.com/nextflow-io/nextflow/pull/7696).
- [nextflow-io/language-server#191](https://github.com/nextflow-io/language-server/pull/191): the same bundle for the language server, whose AOT cache is created on the user's machine.
- [JEP 483](https://openjdk.org/jeps/483), [JEP 514](https://openjdk.org/jeps/514), [JEP 515](https://openjdk.org/jeps/515).
