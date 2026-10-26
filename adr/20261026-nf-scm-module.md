# Extract SCM support into the nf-scm module

- Authors: Jorge Ejarque
- Status: proposed
- Date: 2026-10-26
- Tags: scm, modularization, build, platform

Technical Story: supersedes [#6719](https://github.com/nextflow-io/nextflow/pull/6719)

## Summary

Move the `nextflow.scm` package and `Manifest` into a new `nf-scm` module that depends only on `nf-commons` and JGit. The config parser stays in the `nextflow` runtime, and `nf-scm` finds it at runtime through a `ServiceLoader` hook. Library consumers such as Seqera Platform can then use SCM support without the whole runtime.

## Problem Statement

Seqera Platform depends on `io.nextflow:nextflow`, which is the whole runtime: executors, config and script parsers, GPars, Kryo, LevelDB, grengine, jline and more. Platform uses about 25 classes from it. Almost all of them are SCM code (`AssetManager`, `RepositoryProvider` and its implementations, `ProviderConfig`, `Manifest`, SCM exceptions) plus a few small helpers (`BashTemplateEngine`, `MustacheTemplateEngine`, `ConfigHelper`). The rest of what Platform uses already lives in `nf-commons`.

The `nextflow.scm` package already gets almost all of its imports from `nf-commons`. Only two links point into the runtime:

1. `AssetManager.getScriptFile()` creates a `nextflow.script.ScriptFile`, a runtime class.
2. `AssetManager.getManifest0()` and `ProviderConfig` parse Nextflow config text with `ConfigParserFactory`, which brings in the whole config parser (v1 and v2).

The second link is the hard one. Two files need the parser:

- **The repository `nextflow.config`**, to read the `manifest` block (main script, default branch, submodules, home page).
- **The SCM file** (`~/.nextflow/scm` or `NXF_SCM_FILE`), which is evaluated as Groovy config. It can read environment variables (`env('NAME')`) and contain expressions, so it cannot be read without the evaluating parser.

Platform no longer uses the config parser. It reads the pipeline config without compiling it, and overrides `getManifest0()` so that repository code is never evaluated on the backend. It also passes the SCM configuration as a `Map`, so `ProviderConfig` never parses text there.

## Goals or Decision Drivers

- **Smaller dependency for library consumers**: Platform should get SCM support without the runtime.
- **No config parser in `nf-scm`**: the module must not ship code that evaluates repository config files, since Platform deliberately avoids it.
- **No behaviour change for Nextflow users**: the manifest and the SCM file are parsed exactly as today.
- **No source change for Platform and plugins**: keep the package names (`nextflow.scm`, `nextflow.config.Manifest`), and keep `nextflow` re-exporting the new module.
- **Small change**: mostly class moves, so the change is easy to review and rebase.

## Non-goals

- **Moving SCM support out of the Nextflow repository** or removing the `nf-scm` → `nf-commons` dependency. This would make Platform independent of Nextflow versions, but it needs larger changes on both sides (see Rationale).
- **Reading the manifest without evaluating it in Nextflow.** It would remove the parser need for the manifest, but it changes behaviour, so it is left for a separate change.
- **Migrating Platform.** Platform switches to `nf-scm` in its own change, after a Nextflow release that includes the module.

## Considered Options

The class moves are the same in every option. The options differ in how `nf-scm` gets the config parser.

1. **Service lookup**: the parser stays in `nextflow`, and `nf-scm` finds it with `java.util.ServiceLoader`.
2. **Move the config and secrets packages too**: new `nf-config` and `nf-secrets` modules, and `nf-scm` depends on `nf-config` (the approach of [#6719](https://github.com/nextflow-io/nextflow/pull/6719)).
3. **Move the parsing out of `nf-scm`**: `nf-scm` keeps only the models, and the code that parses files moves to `nextflow`.

## Pros and Cons of the Options

### Option 1: Service lookup

`nf-commons` defines a one-method interface next to the existing `ConfigParser` interface. `nextflow` implements it and registers the implementation in `META-INF/services`.

```java
// nf-commons
public interface ConfigParserProvider {
    ConfigParser create();
    static ConfigParser createParser() { /* ServiceLoader lookup, error if none found */ }
}
```

```groovy
// nextflow
class DefaultConfigParserProvider implements ConfigParserProvider {
    ConfigParser create() { ConfigParserFactory.create() }
}
```

- Good, because `nf-scm` has no compile or runtime dependency on the parser.
- Good, because no caller changes: `AssetManager` and `ProviderConfig` call `ConfigParserProvider.createParser()` where they called `ConfigParserFactory.create()`.
- Good, because it follows an existing pattern: `nf-lang` loads `RemoteModuleResolver` the same way.
- Good, because the fat jar already merges service files (`mergeServiceFiles()` in `nf-cli-v1`).
- Bad, because a missing registration is only found at runtime. A test in `nextflow` checks that the provider is found, and reading a manifest fails with an error instead of returning an empty manifest.

### Option 2: Move the config and secrets packages too

Three new modules: `nf-secrets` (the `secret` package), `nf-config` (the whole `config` package, including `ConfigBuilder` and the v1/v2 parsers) and `nf-scm` (depends on `nf-config`). `ConfigBuilder` is decoupled from the CLI with new options interfaces, and `NF`/`HistoryFile` usages are removed from it.

- Good, because `AssetManager` keeps calling `ConfigParserFactory` directly, with no lookup.
- Good, because `nf-config` could be reused by tools that read Nextflow config without the runtime.
- Bad, because `nf-scm` brings the config parsers, grengine and snakeyaml to every consumer. Platform no longer uses them, and does not want code that evaluates repository config on its classpath.
- Bad, because it is a much larger change (122 files in #6719), touching `ConfigBuilder` and many CLI commands.
- Bad, because it was designed when Platform still used `ConfigBuilder`, which is no longer true.

### Option 3: Move the parsing out of nf-scm

Split `ProviderConfig` into a model (in `nf-scm`) and a factory that reads and parses the SCM file (in `nextflow`). Take the manifest handling out of `AssetManager`: either make `getManifest0()` abstract with a runtime subclass, or pass the manifest in from outside.

- Good, because the dependency direction is explicit, with no runtime lookup.
- Good, because it moves towards an `AssetManager` that does not depend on the config parser at all, which a later move out of the Nextflow repository would need.
- Bad, because the `AssetManager` constructors that read the SCM file are called in 10 places (CLI commands, `K8sDriverLauncher`, `nf-tower`). All of them, and any external plugin calling `new AssetManager(...)`, must change.
- Bad, because `AssetManager` reads the manifest in many places, including its constructor, so taking it out is a redesign rather than a move.

## Solution or decision outcome

Option 1: a new `nf-scm` module, with the config parser found at runtime through `ConfigParserProvider`.

## Rationale & discussion

### Module structure

```
nf-cli-v1 ──► nextflow (runtime) ──► nf-scm ──► nf-commons ──► nf-lang
                                       └──► JGit
```

- **`nf-scm`**: the whole `nextflow.scm` package and `nextflow.config.Manifest`. It has its own `META-INF/extensions.idx` that registers `RepositoryFactory`.
- **Moved to `nf-commons`**, because Platform uses them and they have no runtime dependencies:
  - `BashTemplateEngine` and `MustacheTemplateEngine`;
  - `ConfigHelper` and `ConfigClosurePlaceholder`;
  - the `ConfigParser` interface, and the new `ConfigParserProvider`;
  - the exceptions `ConfigParseException`, `HttpResponseLengthExceedException`, `MissingCredentialsException`, `RateLimitExceededException` and `AmbiguousPipelineNameException`.

  `nf-commons` now also exports `groovy-json` and `snakeyaml`, which `ConfigHelper` needs. Both are already on every Nextflow classpath.
- **`nextflow`** declares `api project(':nf-scm')`, so plugins and the CLI see the same classes as before.

Package names do not change, so Platform and plugins need no source changes.

### The ScriptFile link

`AssetManager.getScriptFile()` is replaced by a static factory in the runtime, `ScriptFile.of(AssetManager, String)`. `CmdRun` and `CmdPull` are its only callers. This is a small API break for external code that calls `getScriptFile()`. No plugin in this repository does.

### The config parser link

The service lookup was chosen because it keeps the parser out of `nf-scm` (the main driver) without changing any caller, and it follows a pattern already used in the codebase.

Option 2 fails the main driver: it ships the parsers to every consumer. Option 3 meets the drivers, but it is a redesign of `AssetManager` with breaking constructor changes, which is out of scope for a module extraction.

Three details of the lookup:

- **Class loader**: it uses `ConfigParserProvider`'s own class loader, not the thread context one, so code running in plugin threads finds the service too.
- **Fail loudly**: `getManifest0()` gets the parser outside its `try` block. Without this, a missing parser would be caught and turned into an empty manifest. Nextflow would then silently use `main.nf`, the remote default branch and no submodules.
- **Where the SCM file fails**: a missing parser fails only when an SCM file exists. When the file is absent, `ProviderConfig` returns an empty configuration without parsing.

A class-name lookup (`Class.forName('nextflow.config.ConfigParserFactory')`) was also tried. It was rejected because a rename would still compile and only fail at runtime, and the dependency would be hidden in a string.

### Tests

The SCM tests move to `nf-scm`. They use `testImplementation(testFixtures(project(':nextflow')))` for `TemporaryPath` and for the runtime parser. This is test scope only, and does not change what `nf-scm` brings to consumers.

These tests stay in `nextflow`:

- `ProviderPathTest`, because it uses `ConfigBuilder`;
- `ScriptFileTest`, for `ScriptFile.of`;
- `DefaultConfigParserProviderTest`, which checks that the provider is found when the runtime is on the classpath.

### Packaging

No packaging change is needed. The distribution takes `nf-cli-v1`'s dependencies, and `nextflow` re-exports `nf-scm`, so the module is included. The fat jar merges `META-INF/services` files and appends all `extensions.idx` files, so both `DefaultConfigParserProvider` and `RepositoryFactory` are registered.

### Impact on consumers

For a consumer that depends on `nf-scm` instead of `nextflow` (measured with Gradle's runtime classpath):

- external libraries drop from 68 to 38;
- among the libraries removed: GPars, Kryo, LevelDB, grengine, jline and npr-client;
- the `nextflow` and `nf-httpfs` jars are no longer included.

A consumer that relied on those libraries arriving through `nextflow` must now declare them. For Platform, that means `groovy-xml` and `nf-httpfs`.

### Future work

- Read the manifest without evaluating it, from the syntax tree built by `nf-lang`. Nextflow would then not run repository code to read it, as Platform already does. The SCM file would still need the evaluating parser.
- Moving SCM support out of the Nextflow repository would let Platform stop following Nextflow versions. It needs Option 3's redesign of `AssetManager`, and removing the `nf-commons` dependency (`Const`, `SysEnv`, `AbortOperationException`, `RetryConfig`, `Plugins`). On the Platform side it needs `PipelineAssets` and `StudioAssetsManager` to stop subclassing `AssetManager`, and `SsrfGuardedRepositoryProvider` to stop relying on package `nextflow.scm` to reach the protected `getAuth()`.

## Links

- Supersedes [#6719](https://github.com/nextflow-io/nextflow/pull/6719): Move shared Platform code to separate modules
- Related to [Separate CLI from Runtime](20260625-separate-cli-from-runtime.md)
