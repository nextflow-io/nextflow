# Standalone bundle with a Java runtime and AOT cache

- Authors: Ben Sherman, Phil Ewels
- Status: draft
- Deciders: TBD
- Date: 2026-10-28
- Tags: launcher, install, java, distribution, performance

Technical Story: [#7696](https://github.com/nextflow-io/nextflow/pull/7696), supersedes [#7697](https://github.com/nextflow-io/nextflow/pull/7697)

## Summary

Each release publishes a standalone bundle per platform, next to the jar, which contains the launcher, a Temurin Java 25 runtime, the jar and a JVM AOT cache trained in CI. With `NXF_JAVA_BUNDLE=true`, the installed `nextflow` script hands off to the bundle of the selected `NXF_VER`, downloading it when needed. Nextflow then starts faster, and no longer depends on the Java installed by the user. Without it, nothing changes.

## Problem Statement

Nextflow runs on whatever Java it finds, which causes two problems:

- **Installation.** Users must install a compatible Java first. A Java upgrade can break a working install: Homebrew had to pin `openjdk@25` after JDK 27 failed with `Unsupported class file major version 71`.
- **Startup.** Most of the startup of a short run is the JVM loading and linking ~6,000 classes. Java 25 can cache this work ahead of time (JEPs 483, 514, 515), but the cache is only valid for the exact runtime, jar and JVM options it was created with, so it can't be shipped for an arbitrary user Java.

`NXF_VER` adds a constraint: the launcher can run any version from 24.07.0-edge. A runtime must match the selected version, not the launcher. Versions 24.10 and 25.04 fail on Java 25 (`Unsupported class file major version 69`, from Groovy's ASM), while 25.10 and later work. A future runtime upgrade will break older versions in the same way.

## Goals or Decision Drivers

- Each version runs with the runtime and JVM options it was built and tested with, including versions selected with `NXF_VER`.
- A launcher installed today can run versions released later, even if they need a newer Java.
- The AOT cache works on first use, including CI jobs and containers with a fresh `NXF_HOME`.
- No change in behaviour for users who don't opt in.
- Nothing gets worse when a bundle is not available (older versions, other platforms, offline): the launcher falls back to the jar and the system Java.

## Non-goals

- Using the bundle by default, or automatically when no Java is installed. Users who want a standalone install will get a separate install path.
- A native binary (GraalVM native image): Nextflow compiles scripts and config files at runtime, and pipelines ship arbitrary Groovy code in `lib/`.
- Windows, Intel macOS and musl Linux bundles (they keep the jar and the system Java).

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
- Bad, because creating the cache in a shared `NXF_HOME` needs coordination between parallel launches, and the cache is silently rejected on hosts with 128 GB of RAM or more when it was created on a smaller host (see below).

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
- Bad, because each version is a larger download (~115-120 MB instead of 49 MB for the jar), and runtimes are not shared between versions.

## Solution or decision outcome

Option 4: a bundle per version and platform (linux-x64, linux-arm64, macos-arm64), with the stock Temurin JRE, published by the release next to the jar. The installed launcher hands off to it only with `NXF_JAVA_BUNDLE=true`.

## Rationale & discussion

### Bundle

```
nextflow-<version>-<os>-<arch>/
  bin/nextflow                  launcher (NXF_VER defaults to the bundled version)
  runtime/                      Temurin 25 JRE
  lib/nextflow-<version>-one.jar
  lib/nextflow.aot              AOT cache
```

- **Runtime.** The latest Temurin 25 JRE, resolved at build time from the Adoptium API and verified with its SHA-256 checksum, so each release ships the latest Java security fixes without a change in the repository. The build resolves the version once for every platform, and fails if a platform has no release of it. `-PjreVersion=<version>` reproduces a build with another version, and the repository variable `NXF_BUNDLE_JRE_VERSION` overrides it for a release, e.g. if the latest is not yet available for every platform. The version is recorded in `runtime/release`.
- **AOT cache.** It is trained by running a small pipeline (`bundle/training`) with the bundle launcher, so that it is created with the same JVM options as when it is used. It can only be created by the target runtime, so each platform is built on a native CI runner.
- **Size.** ~115-120 MB per platform (tar.gz). On linux-x64, the runtime is 62 MB, the jar 44 MB and the AOT cache 14 MB.

### Stock JRE or jlink

The startup gains come from shipping a fixed runtime with a prebuilt AOT cache. `jlink` does not change startup (within noise of a full JDK); what it provides is size:

| runtime (linux-x64) | tar.gz |
|---|---|
| `jlink` image with the modules found by `jdeps` | 43 MB |
| Temurin 25 JRE | 62 MB |
| Temurin 25 JDK | 134 MB |

`jlink` was rejected. Its module list must be maintained, and a missing module only fails at runtime, when a class of it is used: the basic list (`jdeps` + `jdk.zipfs`) broke S3 in nf-amazon (`NoClassDefFoundError: jdk/net/Sockets`). Third-party plugins, JDBC drivers and libraries run in the same JVM, and can need modules that no core code or plugin uses (e.g. `jdk.security.auth`, `java.instrument`, `jdk.jfr`, `jdk.crypto.cryptoki`, `java.xml.crypto`). They can't be tested against the module list, and since the runtime is fixed per version, a missing module could only be fixed with a new release. The stock JRE removes this class of failure for ~19 MB more per download, which trimming the jar can recover (see below). It also needs no `jmods`, which Temurin 25 does not ship.

The size of each module, measured with `jlink`, in case it is reconsidered:

| module | tar.gz | used for |
|---|---|---|
| `java.base` | 18.7 MB | required, mostly the JVM (`libjvm.so`, 28 MB raw) |
| `jdk.localedata` | 9.6 MB | formatting with non-English locales; not referenced by the jar (`--include-locales=en` keeps English only) |
| `java.desktop` | 9.0 MB | required for `java.beans` (Groovy, snakeyaml); its native libs (2.6 MB) are never loaded |
| `jdk.charsets` | 1.0 MB | non-standard charsets |
| `java.net.http` | 0.4 MB | required (lib-httpx, npr-client) |
| `java.security.jgss` | 0.3 MB | Kerberos in jgit |
| `jdk.jdwp.agent` | 0.2 MB | `-remote-debug` |
| others | < 0.1 MB each | |

### Trimming the jar (follow-up)

The Nextflow modules and Groovy are 14 MB of the 44 MB jar. The largest dependencies:

| dependency | size | notes |
|---|---|---|
| Bouncy Castle | 12.4 MB | only used by nf-k8s (embedded since 26.08) to read PEM keys, in `SSLUtils.java` |
| jgit | 3.4 MB | `pull` and `clone` |
| Guava | 3.3 MB | used across nf-commons and nextflow |
| Jackson + json-schema-validator | 3.1 MB | used by one class; overlaps with Gson and groovy-json |
| ANTLR 4 tool + StringTemplate + ANTLR 3 | 1.8 MB | the grammar generator of nf-lang, leaked onto the runtime classpath by the Gradle `antlr` plugin |
| Woodstox | 1.8 MB | registered as the `XMLInputFactory`, which the JDK provides |
| GPars | 1.7 MB | |
| Apache Commons (compress, lang3, io, codec) | 2.7 MB | |
| others (logback, javax.mail, jsoup, kryo, leveldb, ...) | ~5 MB | |

The cheapest reductions are reading PEM keys with the JDK in nf-k8s instead of Bouncy Castle (12.4 MB) and removing the ANTLR tool from the runtime classpath (1.8 MB). Beyond that, each dependency would need to be replaced or removed, which would also shrink the AOT cache. These reductions also apply to the jar of the other distributions.

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

### Opt-in

The bundle is used only with `NXF_JAVA_BUNDLE=true`. Using it by default would silently change the behaviour of existing installs:

- `JAVA_HOME` and the `java` on the `PATH` would be ignored, only `NXF_JAVA_HOME` would select another Java.
- The trust store would change (see below): CAs imported into the user's own JDK would no longer be trusted.
- Java security fixes would no longer come from OS updates: each version runs the runtime it was released with. Old versions, e.g. pinned with `NXF_VER`, keep an old runtime.
- Packages that run `nextflow -download` at build time (bioconda) would ship the bundle of the build platform.

Opt-in leaves existing installs unchanged, and makes these trade-offs a choice. A standalone install path, which installs the launcher with `NXF_JAVA_BUNDLE=true` for users without Java, is a follow-up. The default could be revisited in a stable release, with migration notes, once the bundle has been used in practice.

### Launcher

- When the script is in a bundle that contains the jar for `NXF_VER`, it uses the bundled runtime, jar and AOT cache, unless `NXF_JAVA_HOME` is set.
- Otherwise, with `NXF_JAVA_BUNDLE=true`, on a supported platform and without `NXF_JAVA_HOME`, the installed script hands off to the bundle of `NXF_VER` in `$NXF_HOME/framework/<version>/`, downloading it from `$NXF_BASE/v<version>/` (next to the jar) if it is not installed, even if the jar is. If it is not available (e.g. older versions, or offline), the launcher writes a marker for this version and platform, so that later launches don't retry, and uses the jar and the system Java. Deleting the marker retries. On other platforms, it prints a note and uses the jar.
- The decision is per platform, since an `NXF_HOME` can be shared across platforms (e.g. x86 login nodes and arm64 compute nodes, or a macOS `NXF_HOME` mounted into a Linux container). Without the jar, the launcher uses the jar of a bundle of another platform with the system Java, instead of downloading it, with or without `NXF_JAVA_BUNDLE`.
- The bundle is extracted into a temporary directory, removed if the download is interrupted, and installed with an atomic symlink, so that concurrent first launches install it once. It is only installed if its checksum matches and its runtime runs (it may not, e.g. on a `noexec` mount). The download times out after 10 s without a connection, and reports to stderr, so that the output of commands such as `nextflow config` is not affected.
- The hand-off happens after `-download` and `self-update` are handled by the installed script. The bundle launcher is started with a marker in its environment, which it removes, and fails if it is not in a bundle, to avoid a loop. `NXF_VER` is not exported, and the bundle launcher does not export the `JAVA_CMD` of its runtime, so that tasks (e.g. a child `nextflow` of another version) do not inherit them.
- Mirrors and air-gapped installs can set `NXF_JAVA_BUNDLE_BASE` to another host with the same layout (`<base>/v<version>/<bundle>.tar.gz` and `.sha256`), including a `file://` URL. Mirrors of `NXF_BASE` get the bundles too.
- `NXF_JVM_ARGS=-XX:AOTMode=off` disables the AOT cache.

### TLS trust

The bundled runtime has its own CA certificates, so it would not trust the CAs installed in the OS, e.g. by a corporate proxy that inspects TLS. Unless `NXF_OPTS` or `NXF_JVM_ARGS` set a trust store, the launcher uses the trust store of the OS, as the Java of the OS would:

- **Linux:** the Java trust store of the distro, `/etc/pki/java/cacerts` (RHEL and derivatives, generated by `update-ca-trust`) or `/etc/ssl/certs/java/cacerts` (Debian and Ubuntu, only if `ca-certificates-java` is installed), when it exists. Otherwise the bundled CA certificates are used: on Debian and Ubuntu without `ca-certificates-java`, CAs added with `update-ca-certificates` are only in the PEM bundle, which Java can't use as a trust store.
- **macOS:** the keychain (`-Djavax.net.ssl.trustStoreType=KeychainStore-ROOT`, JDK 23+), which includes the root certificates installed by device management. Not yet tested on a managed Mac.

Users who added CAs to the trust store of their own JDK need `NXF_JAVA_HOME`, or a trust store in `NXF_OPTS`.

### AOT cache validity

The JVM ignores the cache silently (the launcher sets `-Xlog:aot=off`) when:

- the jar mtime differs from training. The bundle archive keeps the file timestamps, which Gradle resets by default for reproducible archives.
- the module options differ, e.g. an extra `--add-opens` in `NXF_OPTS`.
- compressed oops are disabled, which the JVM does for a heap of 32 GB or more, i.e. by default on hosts with 128 GB of RAM or more. The launcher does not force compressed oops, since that would cap the heap below 32 GB (also with `-XX:MaxRAMPercentage`): on these hosts the cache is not used unless the heap is set below 32 GB.
- ZGC (until JDK 26) or compact object headers are used.

G1 or Parallel GC instead of Serial, the C2 compiler, `-Xmx` up to 32 GB, `-D` properties and moving the bundle to another directory are all fine.

### Results

Measured on linux-x64 (16 CPUs), mean of 5 runs, the same jar and launcher options (C1 compiler, serial GC), with the earlier `jlink` runtime (the stock JRE starts the same):

| | full JDK | bundle runtime | bundle runtime + AOT cache |
|---|---|---|---|
| `nextflow info` | 648 ms | 718 ms | 459 ms |
| hello | 1,725 ms | 1,788 ms | 945 ms |
| pipeline with 10 tasks | 1,546 ms | 1,642 ms | 987 ms |
| nf-core/rnaseq `--help` | 8,582 ms | 8,375 ms | 7,553 ms |

The cache saves a roughly constant 0.5-0.8 s per launch, and does not change memory usage or the throughput of long runs. On large pipelines, startup is dominated by script compilation (~5.4 s for rnaseq), which the cache does not cover.

### Release

- The release workflow builds the bundles of every platform on native runners (the `Bundle` workflow, called on `[release]` commits), resolving the JRE version once, and the release job requires them: a release can't be published without its bundles. `release.sh` checks them before anything is published.
- `deploy` uploads them with the other release files to `s3://www2.nextflow.io/releases/v<version>/`, served by the CDN at `www.nextflow.io/releases`. It refuses to deploy a version that exists, so published bundles are never replaced. `NXF_RELEASE_BUNDLES=none` deploys a release without bundles, e.g. by hand.
- The jar in the bundle is built from the same commit and build number as the published jar, in another job, so it is not byte-identical to it.

### Distribution and packaging

- **The `-dist` release.** It is a single file with the launcher and the jar, which needs a system Java. It is still needed for platforms without a bundle, and Homebrew builds it from source with its own `openjdk`. It should be kept. The bundle could replace it for offline installs on supported platforms.
- **bioconda** installs the launcher and runs `nextflow -download` at build time, with `NXF_DIST` set inside the package. It is unaffected, since the bundle is opt-in.
- **Docker image.** It is based on `amazoncorretto:25-al2023` (glibc), so it could use the bundle of its platform instead of installing a JDK.
- **Package managers.** A package manager owns the Nextflow version, the Java dependency and updates. Packages could set a single variable, e.g. `NXF_INSTALLER=conda`, which makes `self-update` print the upgrade command of the package manager instead of overwriting a file that it owns. Not implemented yet.

### Testing

In bare containers (linux-x64), with the bundles served from a local mirror:

| Scenario | Result |
|---|---|
| Rocky 9, CentOS 7 (glibc 2.17), Debian slim, Ubuntu 24.04, no Java installed | bundle installed on first launch, pipeline runs with the AOT cache (`-XX:AOTMode=on`) |
| Alpine (musl) | no bundle, same error as today when no Java is installed |
| Debian slim with only `wget`, or neither `curl` nor `wget`, `file://` mirror | works |
| No network, or dropped packets (blackhole IP) | falls back to the jar immediately, or after the 10 s timeout, and does not retry |
| Read-only `NXF_HOME` with the bundle installed | works |
| `NXF_HOME` on a `noexec` mount | bundle not installed, jar and system Java used |
| Private CA added with `update-ca-trust` (Rocky 9) | HTTPS works through the OS trust store; fails with the bundled CA certificates only |
| 4 concurrent first launches | all succeed, one bundle installed |
| Download interrupted (SIGINT, SIGTERM) | no temporary directory left |
| Only the bundle of another platform installed | its jar is used with the system Java, no download |
| `NXF_VER` without a bundle (25.10.9) | jar and system Java used, no error |

The container scenarios were run with the earlier `jlink` runtime; the launcher behaviour does not depend on it. CI builds each bundle on its platform, runs the training pipeline from an extracted copy with `-XX:AOTMode=on`, reads a file over HTTPS, and tests the hand-off of the installed launcher to the bundle from a `file://` mirror, with a subset of the integration tests on linux-x64.

### Open questions

- The standalone install path for users without Java.
- `self-update` run from an extracted bundle replaces the bundle's launcher with the plain one. It keeps working, through the hand-off, but the bundle is no longer self-contained.
- Old versions accumulate in `$NXF_HOME/framework` as today; `self-update -clean` removes them, including their bundles.
- The training pipeline covers a typical run. Commands such as `lint` load classes that it does not cover, which are loaded as usual.

## Links

- Supersedes [#7697](https://github.com/nextflow-io/nextflow/pull/7697) (launcher-managed runtime) and the local AOT cache of the original [#7696](https://github.com/nextflow-io/nextflow/pull/7696).
- [nextflow-io/language-server#191](https://github.com/nextflow-io/language-server/pull/191): the same bundle for the language server, whose AOT cache is created on the user's machine.
- [JEP 483](https://openjdk.org/jeps/483), [JEP 514](https://openjdk.org/jeps/514), [JEP 515](https://openjdk.org/jeps/515).
