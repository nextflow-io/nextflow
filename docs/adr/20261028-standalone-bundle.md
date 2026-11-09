# Standalone bundle with a Java runtime and AOT cache

- Authors: Ben Sherman, Phil Ewels
- Status: accpeted
- Date: 2026-10-28
- Tags: launcher, install, java, distribution, performance

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

- Using the bundle by default, or automatically when no Java is installed. Users who want it install Nextflow with `NXF_JAVA_BUNDLE=true` (see below).
- A native binary (GraalVM native image): Nextflow compiles scripts and config files at runtime, and pipelines ship arbitrary Groovy code in `lib/`.
- Windows, Intel macOS and musl Linux bundles (they keep the jar and the system Java).

## Considered Options

The bundle answers two questions:

1. **Runtime:** the stock Temurin JRE, or a `jlink` image with only the modules that Nextflow uses.
2. **AOT cache:** created in CI and shipped with each version, or created on the user's machine.

## Pros and Cons of the Options

### Stock JRE

- Good, because every JDK module is available: third-party plugins, JDBC drivers and libraries run in the same JVM, and can need modules that no core code or plugin uses (e.g. `jdk.security.auth`, `java.instrument`, `jdk.jfr`, `jdk.crypto.cryptoki`, `java.xml.crypto`).
- Good, because there is no module list to maintain, and it needs no `jmods`, which Temurin 25 does not ship.
- Bad, because it is 19 MB larger than a `jlink` image, which trimming the jar can recover (see below).

### `jlink` image

| runtime (linux-x64) | tar.gz |
|---|---|
| `jlink` image with the modules found by `jdeps` | 43 MB |
| Temurin 25 JRE | 62 MB |
| Temurin 25 JDK | 134 MB |

- Good, because it is smaller. Startup is the same (within noise of a full JDK): the gains come from the AOT cache.
- Bad, because a missing module only fails at runtime, when a class of it is used: the basic list (`jdeps` + `jdk.zipfs`) broke S3 in nf-amazon (`NoClassDefFoundError: jdk/net/Sockets`).
- Bad, because third-party code can't be tested against the module list, and since the runtime is fixed per version, a missing module could only be fixed with a new release.

### AOT cache created in CI

The cache is only valid for the exact runtime build, jar and JVM options it was created with, so the runtime must ship with each version: a bundle per version and platform, to which the installed launcher hands off.

- Good, because the cache works on first use, including CI jobs and containers with a fresh `NXF_HOME`, and each install runs the artifact that was tested.
- Good, because each version runs with the runtime and JVM options it was built with, whatever launcher is installed: a launcher installed today can run a later version that needs a newer Java.
- Good, because the bundle is self-contained: it can be extracted and run directly, e.g. in container images or air-gapped sites.
- Bad, because runtimes are not shared between versions: each version is 110 MB instead of 38 MB for the jar.
- Bad, because the installed script must keep supporting the bundle URL and layout for as long as old launchers exist, and the jar and system Java for versions and platforms without a bundle.

### AOT cache created on the user's machine

The launcher downloads a runtime for `NXF_VER`, shared between versions that use the same Java. The first run records a training profile, and a later launch creates the cache.

- Good, because runtimes are shared between versions, and downloaded separately from the jar.
- Bad, because the first runs don't benefit, and CI jobs and containers with a fresh `NXF_HOME` never do.
- Bad, because creating the cache in a shared `NXF_HOME` needs coordination between parallel launches, and a cache created on a smaller host is silently rejected on hosts with 128 GB of RAM or more (see below).
- Bad, because the launcher must choose the runtime of each version: an installed launcher only knows the versions released before it, and the JVM options come from the installed launcher rather than from the version that runs.

## Solution or decision outcome

The stock Temurin JRE, and an AOT cache created in CI: each release publishes a bundle per platform (linux-x64, linux-arm64, macos-arm64) next to the jar. The installed launcher hands off to it only with `NXF_JAVA_BUNDLE=true`.

## Rationale & discussion

### Bundle

```
nextflow-<version>-<os>-<arch>/
  bin/nextflow                  launcher (NXF_VER defaults to the bundled version)
  runtime/                      Temurin 25 JRE
  lib/nextflow-<version>-one.jar
  lib/nextflow.aot              AOT cache
```

- **Runtime.** The latest Temurin 25 JRE, resolved at build time from the Adoptium API and verified with its SHA-256 checksum, so each release ships the latest Java security fixes without a change in the repository. The version is resolved once for every platform: if the latest release differs between them, the newest release that every platform has is used. `-PjreVersion=<version>` or the repository variable `NXF_BUNDLE_JRE_VERSION` override it. The version is recorded in `runtime/release`.
- **AOT cache.** It is trained by running a small pipeline (`bundle/training`) with the bundle launcher, so that it is created with the same JVM options as when it is used. It can only be created by the target runtime, so each platform is built on a native CI runner.
- **Size.** 110 MB per platform (tar.gz). On linux-x64, the runtime is 62 MB, the jar 34 MB and the AOT cache 14 MB.

### Trimming the jar (follow-up)

The Nextflow modules and Groovy are 14 MB of the 38 MB jar (raw, Groovy 5). The largest dependencies:

| dependency | size | notes |
|---|---|---|
| jgit | 3.4 MB | `pull` and `clone` |
| Guava | 3.3 MB | used across nf-commons and nextflow |
| Jackson + json-schema-validator | 3.1 MB | used by one class; overlaps with Gson and groovy-json |
| ANTLR 4 tool + StringTemplate + ANTLR 3 | 1.8 MB | the grammar generator of nf-lang, leaked onto the runtime classpath by the Gradle `antlr` plugin |
| Woodstox | 1.8 MB | the `XMLInputFactory` for nf-amazon (#6369): the JDK one retains parsed responses when the AWS SDK creates one per virtual thread (> 500 MB for 15K uploads) |
| GPars | 1.7 MB | |
| Apache Commons (compress, lang3, io, codec) | 2.7 MB | |
| others (logback, javax.mail, jsoup, kryo, leveldb, ...) | ~5 MB | |

The possible reductions, which also apply to the jar of the other distributions:

- **A hand-written parser for nf-lang** (under consideration) removes ANTLR entirely. Otherwise, excluding the `antlr` configuration from `api` drops the tool and keeps the runtime.
- **Woodstox** could move to nf-amazon, if the plugin sets the thread context classloader around SDK calls (`XMLInputFactory.newInstance()` looks it up there). Dropping it would bring back the heap growth.
- **javax.mail and jsoup** (1.3 MB) are only used by notifications, which could move to a plugin.
- **json-schema-validator** saves only 0.6 MB: Jackson is also needed by npr-api, lib-httpx and groovy-yaml.

Beyond these, each dependency would need to be replaced, which would also shrink the AOT cache.

### Trimming startup latency (follow-up)

The AOT cache removes most of the class loading, but startup is still most of a short run. Here, startup is the time to the first task, i.e. from launching `nextflow` to the first task being submitted.

Breakdown on master (d2b807a07, after the Groovy 5 upgrade) with the linux-x64 bundle, median of 20 runs on 16 CPUs. hello is a single process that runs `echo`, and train is the training pipeline of the AOT cache:

| phase | hello | train |
|---|---:|---:|
| bash launcher | 27 ms | 27 ms |
| JVM creation (map the AOT cache) | 66 ms | 68 ms |
| `main` to the first log line: Groovy runtime, `Launcher`, CLI parsing, logging | 324 ms | 329 ms |
| `run` command setup | 44 ms | 43 ms |
| plugin manager, config load and parse | 34 ms | 213 ms |
| secrets, `Session` construction | 86 ms | 32 ms |
| `Session` init: history, cache, work dir, executors, observers | 86 ms | 78 ms |
| script parse and compile | 128 ms | 128 ms |
| script run up to the first task: process and workflow setup | 92 ms | 213 ms |
| first task: create, hash, build the wrapper, submit | 74 ms | 86 ms |
| **time to first task** | **981 ms** (of 1,188 ms) | **1,218 ms** (of 1,519 ms) |
| of which GC pauses | 47 ms | 114 ms |

The phases don't add up to the total because each is a median. The phases run one after another on the main thread, so each one adds to the time to first task. These numbers don't compare directly with the Results above, which were measured with 26.10.3 (Groovy 4).

### Removing the `-dist` release

With the bundle, Nextflow has two standalone releases: `-dist`, the launcher with the jar appended, and the bundle. `-dist` is marked as legacy with this change and will be removed by 27.04.

The bundle does not cover every use of `-dist`:

| | `-dist` | bundle |
|---|---|---|
| platforms | any with Java 17+, including musl (Alpine), macOS x64, ppc64le | linux-x64, linux-arm64 (glibc), macos-arm64 |
| Java | the system Java, which some sites and Homebrew require | its own |
| size | 38 MB | 110 MB |

What only `-dist` provides is a single file that runs offline and is published on GitHub releases (the jar is only on the CDN). Its replacement is the bundle, or the launcher and the jar for users who need the system Java.

Until 27.04:

- Publish `nextflow-<ver>-one.jar` on GitHub releases, next to `-dist`.
- Mark `-dist` as legacy in the docs and release notes, and document the offline install with the launcher and the jar in `$NXF_HOME/framework/<ver>/`, or the bundle.
- **Homebrew** builds `-dist` from source (`gradle pack`) and installs it as `libexec/nextflow`: its formula must install the launcher and the jar instead.
- **nf-core/setup-nextflow** downloads `nextflow-<ver>-dist` from GitHub for 24.10 and later when `all: true`: it must download the jar instead.

In 27.04, delete `packDist` in `packing.gradle` and the `NXF_PACK=dist` branches of the launcher (`NXF_PACK` is not read outside the launcher), and replace the "Standalone distribution" section of the docs with a migration note. User scripts and air-gapped installs that still download `-dist` break then.

### Opt-in

The bundle is used only with `NXF_JAVA_BUNDLE=true`. Using it by default would silently change the behaviour of existing installs:

- `JAVA_HOME` and the `java` on the `PATH` would be ignored, only `NXF_JAVA_HOME` would select another Java.
- The trust store would change (see below): CAs imported into the user's own JDK would no longer be trusted.
- Java security fixes would no longer come from OS updates: each version runs the runtime it was released with. Old versions, e.g. pinned with `NXF_VER`, keep an old runtime.
- Packages that run `nextflow -download` at build time (bioconda) would ship the bundle of the build platform.

Opt-in leaves existing installs unchanged, and makes these trade-offs a choice. The default could be revisited in a stable release, with migration notes, once the bundle has been used in practice.

Users without Java install Nextflow with `curl -s https://get.nextflow.io | NXF_JAVA_BUNDLE=true bash`. The installer then writes `true` as the default of `NXF_JAVA_BUNDLE` into the installed launcher, so that the setting does not depend on the environment, and `self-update` keeps it. The launcher of a bundle has the same default, so that an extracted bundle also runs other versions (`NXF_VER`) with their bundle. `NXF_JAVA_BUNDLE=false` still selects the system Java.

### Launcher

- When the script is in a bundle that contains the jar for `NXF_VER`, it uses the bundled runtime, jar and AOT cache, unless `NXF_JAVA_HOME` is set.
- Otherwise, with `NXF_JAVA_BUNDLE=true`, on a supported platform and without `NXF_JAVA_HOME`, the installed script hands off to the bundle of `NXF_VER` in `$NXF_HOME/framework/<version>/`, downloading it from `$NXF_BASE/v<version>/` (next to the jar) if it is not installed, even if the jar is. With `NXF_OFFLINE=true`, it is not downloaded. On other platforms, the launcher prints a note and uses the jar.
- If the bundle can't be used (not on the server, i.e. HTTP 403 or 404 from the CDN, e.g. for older versions; a bad checksum; or a runtime that doesn't run), the launcher writes a marker for this version and platform, so that later launches don't retry, and uses the jar and the system Java. Deleting the marker retries. A failed download (e.g. no network) is not marked and is retried at the next launch, so that an install without Java recovers.
- The decision is per platform, since an `NXF_HOME` can be shared across platforms (e.g. x86 login nodes and arm64 compute nodes, or a macOS `NXF_HOME` mounted into a Linux container).
- When the launcher does not use a bundle and the jar of `NXF_VER` is not installed, it uses the jar of an installed bundle with the system Java, instead of downloading the jar (which ignores `NXF_OFFLINE`). This happens e.g. when `NXF_JAVA_HOME` is set after an install without Java, or in an `NXF_HOME` shared with a platform that has no bundle.
- The bundle is installed with an atomic symlink, so that concurrent first launches install it once, and only if its checksum matches and its runtime runs (it may not, e.g. on a `noexec` mount).
- The bundle launcher is started with a marker in its environment, and fails if it is not in a bundle, to avoid a loop. Neither launcher exports `NXF_VER` or the bundled `JAVA_CMD`, so that tasks (e.g. a child `nextflow` of another version) do not inherit them.
- `self-update` fails in a bundle, e.g. an extracted one: replacing its launcher would leave its runtime and jar unused, and install the new version in `$NXF_HOME` instead. An extracted bundle is updated by extracting the bundle of the new version.
- Old versions accumulate in `$NXF_HOME/framework` as today; `self-update -clean` removes them, including their bundles.
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

The training pipeline covers a typical run: commands such as `lint` load classes that it does not cover, which are loaded as usual.

G1 or Parallel GC instead of Serial, the C2 compiler, `-Xmx` up to 32 GB, `-D` properties and moving the bundle to another directory are all fine.

### Results

Measured on linux-x64 (16 CPUs) with the linux-x64 bundle of 26.10.3, median of 10 runs. The full JDK is Temurin 25.0.5, the same build as the bundle runtime, run through the bundle launcher with `NXF_JAVA_HOME`, so the jar and JVM options (C1 compiler, serial GC) are the same:

| | full JDK | bundle runtime | bundle runtime + AOT cache |
|---|---|---|---|
| `nextflow info` | 582 ms | 714 ms | 453 ms |
| hello | 1,386 ms | 1,461 ms | 891 ms |
| pipeline with 10 tasks | 1,445 ms | 1,587 ms | 1,019 ms |
| nf-core/rnaseq `--help` | 9,217 ms | 8,751 ms | 8,823 ms |

The cache saves 0.5-0.6 s on short runs (0.1-0.3 s for `info`), and lowers their peak RSS by ~10% (e.g. 239 MB instead of 264 MB for hello). On large pipelines, startup is dominated by script compilation, which the cache does not cover: on rnaseq, the gain is within noise.

### Release

- The release workflow builds the bundles of every platform on native runners (the `Bundle` workflow, called on `[release]` commits), and the release job requires them: a release can't be published without its bundles, since an install without Java can't run a version without them, and they can't be added afterwards. `release.sh` checks them before anything is published.
- `deploy` uploads them with the other release files to `s3://www2.nextflow.io/releases/v<version>/`, served by the CDN at `www.nextflow.io/releases`. It refuses to deploy a version that exists, so published bundles are never replaced. `NXF_RELEASE_BUNDLES=none` deploys a release without bundles, e.g. by hand.
- The jar in the bundle is built from the same commit and build number as the published jar, in another job, so it is not byte-identical to it.

### Testing

CI builds each bundle on its platform, reads a file over HTTPS, and tests the hand-off of the installed launcher to the bundle from a `file://` mirror with `-XX:AOTMode=on`, with a subset of the integration tests on linux-x64. The failure modes of the launcher (no network, missing bundle, bad checksum, `noexec`, concurrent launches, other distros) were tested by hand in bare containers.

## Links

- [#7696](https://github.com/nextflow-io/nextflow/pull/7696): the implementation.
- [nextflow-io/language-server#191](https://github.com/nextflow-io/language-server/pull/191): the same bundle for the language server, whose AOT cache is created on the user's machine.
- [JEP 483](https://openjdk.org/jeps/483), [JEP 514](https://openjdk.org/jeps/514), [JEP 515](https://openjdk.org/jeps/515).
