# Extract SCM support into the nf-scm module

- Authors: Jorge Ejarque
- Status: proposed
- Date: 2026-10-26
- Tags: scm, modularization, build, platform

Technical Story: supersedes [#6719](https://github.com/nextflow-io/nextflow/pull/6719)

## Summary

Move the `nextflow.scm` package and `Manifest` into a new `nf-scm` module that depends only on `nf-commons` and JGit. `nf-scm` no longer uses the evaluating config parser: the SCM file is parsed in the `nextflow` runtime and passed in as a `Map`, and the manifest is read from the config syntax tree without evaluating it. Library consumers such as Seqera Platform can then use SCM support without the whole runtime.

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
- **No new extension mechanism**: avoid adding a public lookup API whose only purpose is to reach back into the runtime.
- **Keep user-facing behaviour**: the SCM file is parsed exactly as today. The manifest may change only where it holds values that are not literals.
- **No source change for Platform**: keep the package names (`nextflow.scm`, `nextflow.config.Manifest`) and the `AssetManager(String, Map)` constructor, and keep `nextflow` re-exporting the new module.

## Non-goals

- **Moving SCM support out of the Nextflow repository** or removing the `nf-scm` → `nf-commons` dependency. This would make Platform independent of Nextflow versions, but it needs larger changes on both sides.
- **Migrating Platform.** Platform switches to `nf-scm` in its own change, after a Nextflow release that includes the module.

## Considered Options

The class moves are the same in every option. The options differ in how `nf-scm` handles config parsing.

1. **Service lookup**: the parser stays in `nextflow`, and `nf-scm` finds it with `java.util.ServiceLoader`.
2. **Move the config and secrets packages too**: new `nf-config` and `nf-secrets` modules, and `nf-scm` depends on `nf-config` (the approach of [#6719](https://github.com/nextflow-io/nextflow/pull/6719)).
3. **Remove config evaluation from `nf-scm`**: the SCM file is parsed in `nextflow` and passed in as a `Map`, and the manifest is read without evaluating it.

## Pros and Cons of the Options

### Option 1: Service lookup

`nf-commons` defines a `ConfigParserProvider` interface with a static `createParser()` lookup. `ConfigParserFactory` in `nextflow` implements it and is registered in `META-INF/services`, and every caller gets its parser through `ConfigParserProvider.createParser()`.

- Good, because no caller of `AssetManager` changes, and user-facing behaviour does not change at all.
- Good, because it follows an existing pattern: `nf-lang` loads `RemoteModuleResolver` the same way.
- Bad, because it adds a public extension point to `nf-commons` whose only purpose is to reach back into the runtime.
- Bad, because a missing registration is only found at runtime, so it needs a guard test and an explicit failure in `getManifest0()`.
- Bad, because `nf-scm` still evaluates repository config to read a manifest whenever it runs inside Nextflow.

### Option 2: Move the config and secrets packages too

Three new modules: `nf-secrets` (the `secret` package), `nf-config` (the whole `config` package, including `ConfigBuilder` and the v1/v2 parsers) and `nf-scm` (depends on `nf-config`). `ConfigBuilder` is decoupled from the CLI with new options interfaces, and `NF`/`HistoryFile` usages are removed from it.

- Good, because `AssetManager` keeps calling `ConfigParserFactory` directly, with no lookup.
- Good, because `nf-config` could be reused by tools that read Nextflow config without the runtime.
- Bad, because `nf-scm` brings the config parsers, grengine and snakeyaml to every consumer. Platform no longer uses them, and does not want code that evaluates repository config on its classpath.
- Bad, because it is a much larger change (122 files in #6719), touching `ConfigBuilder` and many CLI commands.
- Bad, because it was designed when Platform still used `ConfigBuilder`, which is no longer true.

### Option 3: Remove config evaluation from nf-scm

- **SCM file**: split `ProviderConfig` into the model (stays in `nf-scm`) and a `ProviderConfigParser` (in `nextflow`) that finds, reads and parses the file. Callers parse the file and pass the result to `AssetManager` as a `Map`, as Platform already does.
- **Manifest**: a `ManifestParser` in `nf-scm` reads the `manifest` scope from the syntax tree built by `nf-lang`'s `ConfigParser`, without evaluating it.

- Good, because `nf-scm` contains no code that evaluates config, and needs no lookup or registration.
- Good, because the config parser classes stay in `nextflow`, unchanged.
- Good, because Nextflow no longer runs repository code just to read a manifest (`pull`, `info`, `clone`). Platform already reads the pipeline config this way.
- Good, because the CLI and Platform use `AssetManager` the same way: the caller supplies the SCM configuration.
- Bad, because two `AssetManager` constructors are removed, so 10 calls (CLI commands, `K8sDriverLauncher`, `nf-tower`) and any external code calling them must change.
- Bad, because manifest values that are not literals are no longer read (see below).

## Solution or decision outcome

Option 3: a new `nf-scm` module with no config evaluation. The SCM file is parsed by `ProviderConfigParser` in the runtime, and the manifest is read by `ManifestParser` without evaluating it.

## Rationale & discussion

### Module structure

```
nf-cli-v1 ──► nextflow (runtime) ──► nf-scm ──► nf-commons ──► nf-lang
                                       └──► JGit
```

- **`nf-scm`**: the whole `nextflow.scm` package, `nextflow.config.Manifest` and the new `nextflow.config.ManifestParser`. It has its own `META-INF/extensions.idx` that registers `RepositoryFactory`.
- **Moved to `nf-commons`**, because Platform uses them and they have no runtime dependencies:
  - `BashTemplateEngine` and `MustacheTemplateEngine`;
  - `ConfigHelper` and `ConfigClosurePlaceholder`;
  - the exceptions `HttpResponseLengthExceedException`, `MissingCredentialsException`, `RateLimitExceededException` and `AmbiguousPipelineNameException`.

  `nf-commons` now also exports `groovy-json` and `snakeyaml`, which `ConfigHelper` needs. Both are already on every Nextflow classpath.
- **`nextflow`** declares `api project(':nf-scm')`, so plugins and the CLI see the same classes as before. The new `nextflow.config.ProviderConfigParser` lives here.

Package names do not change, so Platform and plugins need no source changes beyond the `AssetManager` constructors described below.

### The ScriptFile link

`AssetManager.getScriptFile()` is replaced by a static factory in the runtime, `ScriptFile.of(AssetManager, String)`. `CmdRun` and `CmdPull` are its only callers. This is a small API break for external code that calls `getScriptFile()`. No plugin in this repository does.

### The SCM file

`ProviderConfig` keeps the model: the provider attributes, `createFromMap()`, `createDefault()` and `getScmConfigPath()`. Reading and parsing the file (`getDefault()`, `getFromFile()`, `parse()`) moves unchanged to `ProviderConfigParser` in `nextflow`. Its behaviour does not change: the same evaluating parser, the same `env('NAME')` support, and the same handling of a missing file.

`AssetManager` no longer reads the SCM file itself:

- **Removed**: `AssetManager(String, HubOptions)` and `AssetManager(String, String revision, String mainScript, HubOptions)`. This also removes the one-argument form `new AssetManager(name)`.
- **Added**: `AssetManager(String, Map scmConfig, HubOptions, String revision = null, String mainScript = null)`.
- **Unchanged**: `AssetManager(String, Map)`, which Platform uses.

Callers in the runtime pass `ProviderConfigParser.getDefault()`:

```groovy
new AssetManager(pipelineName, ProviderConfigParser.getDefault(), hubOpts, revision, mainScript)
```

### The manifest

`ManifestParser` parses the config text with `nf-lang`'s `ConfigParser`, which only builds the syntax tree, and collects the `manifest.*` assignments, written either as a `manifest { }` block or as dotted names. It keeps literal values: strings, numbers, booleans, lists and maps. Assignments inside other scopes, such as `profiles`, are ignored. Includes are ignored, as before.

The v2 syntax is now the default, and the v2 parser is built on the same `nf-lang` grammar. So any `nextflow.config` that Nextflow can parse is also readable here. The differences are in evaluation only:

| Case | Before | After |
|---|---|---|
| Literal values (almost all manifests) | Read | Read, unchanged |
| Value that is not a literal, e.g. `version = "1.0-${suffix}"` | Evaluated | Ignored (debug log). For `mainScript` or `defaultBranch`, the default is used |
| `nextflow.config` with v1-only syntax (`NXF_SYNTAX_PARSER=v1`) | Read | Syntax error, so an empty manifest and a warning |
| Repository code run to read the manifest | Yes | No |

`getManifest0()` stays `protected`, so existing overrides keep working.

### Tests

The SCM tests move to `nf-scm`, together with the new `ManifestParserTest`. They use `testImplementation(testFixtures(project(':nextflow')))` for the `TemporaryPath` helper only. This is test scope only, and does not change what `nf-scm` brings to consumers.

These tests stay in `nextflow`:

- `ProviderPathTest`, because it uses `ConfigBuilder`;
- `ScriptFileTest`, for `ScriptFile.of`;
- `ProviderConfigParserTest`, for parsing the SCM file.

### Packaging

No packaging change is needed. The distribution takes `nf-cli-v1`'s dependencies, and `nextflow` re-exports `nf-scm`, so the module is included. The fat jar appends all `extensions.idx` files, so `RepositoryFactory` is registered.

### Impact on consumers

For a consumer that depends on `nf-scm` instead of `nextflow` (measured with Gradle's runtime classpath):

- external libraries drop from 68 to 38;
- among the libraries removed: GPars, Kryo, LevelDB, grengine, jline and npr-client;
- the `nextflow` and `nf-httpfs` jars are no longer included.

A consumer that relied on those libraries arriving through `nextflow` must now declare them. For Platform, that means `groovy-xml` and `nf-httpfs`.

## Links

- Supersedes [#6719](https://github.com/nextflow-io/nextflow/pull/6719): Move shared Platform code to separate modules
- Related to [Separate CLI from Runtime](20260625-separate-cli-from-runtime.md)
