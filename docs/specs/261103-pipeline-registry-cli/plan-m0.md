# Pipeline Registry CLI — M0 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Include registry pipelines with `include { workflow as X } from 'scope/name'`, auto-install them from the Git commit pinned by the registry release, and add `nextflow pipeline install|list|remove|search|view`.

**Architecture:**
- **nf-lang** decides whether an include is a pipeline include, and finds its base directory (the nearest pipeline root or vendored module directory). It then calls a new SPI method, `RemoteModuleResolver.resolvePipeline`.
- **`modules/nextflow`** implements that method with a `PipelineManager`. The manager reads the registry through `npr-client`, fetches the pinned commit with JGit, and writes `.pipeline-info`.
- **`nf-cli-v1`** adds a `CmdPipeline` command group, modelled on `CmdModule`.

**Tech Stack:** Groovy/Java 21, Spock, JGit (already a dependency of `modules/nextflow`), `io.seqera:npr-client` (pipeline methods from seqeralabs/nextflow-registry#439).

**Spec:** `docs/specs/261103-pipeline-registry-cli/spec.md`, section M0. **ADR:** `docs/adr/20261103-pipeline-registry-cli.md`.

## Global Constraints

- Pipeline include: the source is `scope/name` and some entry name is `workflow` or `params`. Every other non-local include stays a module include.
- Install location: `<base>/pipelines/<scope>/<name>/`. Info file: `.pipeline-info` (Java `Properties`, like `.module-info`).
- `.pipeline-info` keys: `name`, `version`, `workflowPath`, `source` (`git`), `repositoryUrl`, `revision`, `commitId`, `checksum`, `registryUrl`.
- Checksum: SHA-256 hex over relative path + contents of each file, in sorted order (same scheme as `io.seqera.npr.client.ModuleChecksum.computeDirectory`). It leaves out `.pipeline-info`, `.git/`, and any subdirectory that holds its own `.module-info` or `.pipeline-info`.
- Auto-install picks the latest published release and logs `Installed pipeline <scope/name>@<version>`.
- Never overwrite a locally modified pipeline without `-force`.
- The language server must not install anything. `FallbackRemoteModuleResolver` only looks on disk.
- New source files carry the Apache 2.0 header. Commits are signed off (`git commit -s`). Don't touch `plugins/*/VERSION` or `changelog.txt`.
- The root `build.gradle` `mavenLocal()` line is local only: **never commit it**. Commit the `nprVersion` bump only after npr-client is published to S3.

## Review Focus

1. **A fetch fails halfway** (network error, wrong commit): nothing may be left under `pipelines/`, and a rerun must retry the install. Test: Task 4, `failed install leaves no directory`.
2. **The tag moved since publish** (HEAD of `gitRevision` ≠ `commitId`): install must fail with a clear error. It must not install other code than the release pinned. Test: Task 3, `fails when the revision no longer points to the pinned commit`.
3. **Auto-installing a module inside an installed pipeline** must not turn the pipeline into `modified`. Test: Task 2, `ignores nested dirs with their own info file`.
4. **Names with path tricks** (`../x`, `nf-core/a/b`, upper case): reject them before touching the disk. Test: Task 4, `rejects invalid pipeline names`.
5. **A pipeline include without network**, when the pipeline is already installed: no registry call at all. Test: Task 4, `uses installed pipeline without calling the registry`.

---

### Task 0: Depend on npr-client with pipeline methods

**Files:**
- Modify: `gradle.properties` (`nprVersion`)
- Local only, do not commit: root `build.gradle`, add `mavenLocal()` to the `allprojects { repositories {} }` block

- [ ] **Step 1:** In seqeralabs/nextflow-registry branch `261103-pipeline-client` (#439), run `./gradlew -Psnapshot=true :npr-api:publishToMavenLocal :npr-client:publishToMavenLocal`.
- [ ] **Step 2:** Set `nprVersion=0.29.1-SNAPSHOT` and add the local `mavenLocal()` line.
- [ ] **Step 3:** Run `make compile`, then `./gradlew :nextflow:test --tests 'nextflow.module.*' :nf-cli-v1:test --tests 'nextflow.cli.module.*'`. Expected: BUILD SUCCESSFUL, 0 failures.
- [ ] **Step 4:** No commit yet. The version bump is committed in the last task, once the release exists.

---

### Task 1: Pipeline include resolution in nf-lang

**Files:**
- Modify: `modules/nf-lang/src/main/java/nextflow/module/spi/RemoteModuleResolver.java`
- Modify: `modules/nf-lang/src/main/java/nextflow/module/spi/FallbackRemoteModuleResolver.java`
- Modify: `modules/nf-lang/src/main/java/nextflow/script/control/ModuleResolver.java` (`getIncludeUri`, `resolveInclude`, `visitDependencies`)
- Modify: `modules/nf-lang/src/main/java/nextflow/script/control/ResolveIncludeVisitor.java:99`
- Test: `modules/nf-lang/src/test/groovy/nextflow/module/spi/RemoteModuleResolverTest.groovy`, `modules/nf-lang/src/test/groovy/nextflow/script/control/ResolveIncludeTest.groovy`

**Interfaces:**
- Produces, on `RemoteModuleResolver`:
  - `String PIPELINES_DIR = "pipelines"`
  - `String PIPELINE_INFO_FILE = ".pipeline-info"`
  - `static boolean isPipelineInclude(Collection<String> entryNames)`: true if it contains `"workflow"` or `"params"`
  - `default Path resolvePipeline(String name, Path baseDir)`: returns the path of the pipeline's main script
  - `static Path resolveBaseDir(Path includingFile, Path projectDir)`: new rule (below)
- Produces, on `ModuleResolver`: `static URI getIncludeUri(URI uri, IncludeNode node, Path projectDir)`. It replaces the `String source` overload at all three call sites. Keep the `String` overload, delegating with "not a pipeline", only if something else still calls it.

- [ ] **Step 1: Write the failing tests** in `RemoteModuleResolverTest`:

```groovy
def 'should resolve against the nearest pipeline root'() {
    given:
    def root = projectDir.resolve('pipelines/nf-core/rnaseq')
    Files.createDirectories(root.resolve('subworkflows/local'))
    Files.createFile(root.resolve('.pipeline-info'))
    expect:
    RemoteModuleResolver.resolveBaseDir(root.resolve('main.nf'), projectDir) == root
    RemoteModuleResolver.resolveBaseDir(root.resolve('subworkflows/local/align.nf'), projectDir) == root
}

def 'should prefer a vendored module dir inside a pipeline root'() {
    given:
    def root = projectDir.resolve('pipelines/nf-core/rnaseq')
    def wf = root.resolve('modules/nf-core/bam_sort')
    Files.createDirectories(wf)
    Files.createFile(root.resolve('.pipeline-info'))
    expect:
    RemoteModuleResolver.resolveBaseDir(wf.resolve('main.nf'), projectDir) == wf
}

def 'should resolve a nested pipeline against its own root'() {
    given:
    def outer = projectDir.resolve('pipelines/nf-core/a')
    def inner = outer.resolve('pipelines/nf-core/b')
    Files.createDirectories(inner)
    Files.createFile(outer.resolve('.pipeline-info'))
    Files.createFile(inner.resolve('.pipeline-info'))
    expect:
    RemoteModuleResolver.resolveBaseDir(inner.resolve('main.nf'), projectDir) == inner
}

def 'should detect a pipeline include by its entry names'() {
    expect:
    RemoteModuleResolver.isPipelineInclude(['workflow'])
    RemoteModuleResolver.isPipelineInclude(['params'])
    !RemoteModuleResolver.isPipelineInclude(['FASTQC'])
}
```

In `ResolveIncludeTest`, using the fallback resolver and a temp project:
- `include { workflow as RNASEQ } from 'nf-core/rnaseq'`, with `pipelines/nf-core/rnaseq/main.nf` present (typed, with an entry workflow). Expect no errors, and `RNASEQ` resolves to the entry workflow.
- The same include with no `pipelines/` folder. Expect an error message that contains `Pipeline 'nf-core/rnaseq' not found in 'pipelines' directory`.
- `include { FASTQC } from 'nf-core/fastqc'`, with only `pipelines/nf-core/fastqc/` present. Expect the existing module "not found in 'modules' directory" error. This proves that module includes do not change.

- [ ] **Step 2:** Run `make test module=nf-lang class=RemoteModuleResolverTest` and `make test module=nf-lang class=ResolveIncludeTest`. Expected: the new tests FAIL (missing methods).

- [ ] **Step 3: Implement.**
  - `resolveBaseDir(includingFile, projectDir)`:
    1. `parent` = the including file's directory. Keep the existing null handling.
    2. `root` = the nearest directory, from `parent` upward and stopping at `projectDir` (inclusive), that contains `PIPELINE_INFO_FILE`. If there is none, `root` = `projectDir`.
    3. Return `parent` if `isVendoredModuleDir(parent, root)`, otherwise `root`. `isVendoredModuleDir` already takes a root argument; pass `root` in place of `projectDir`.
  - `resolvePipeline` default and `FallbackRemoteModuleResolver`:
    - Resolve `baseDir/pipelines/<name>`.
    - Reject paths that escape `pipelines/`, as the module code does.
    - If the directory is missing, throw `IllegalStateException("Pipeline '<name>' not found in 'pipelines' directory -- use 'nextflow pipeline install' to install it first")`.
    - Return `<dir>/main.nf`.
  - `getIncludeUri(uri, node, projectDir)`:
    - Local sources: unchanged.
    - Otherwise compute `base`, then call `resolvePipeline(source, base)` when `isPipelineInclude(entry names)`, else `resolve(source, base)`.

- [ ] **Step 4:** Run `make test module=nf-lang`. Expected: all PASS, including the existing `RemoteModuleResolverTest` cases.

- [ ] **Step 5: Commit**

```bash
git add modules/nf-lang
git commit -s -m "Resolve pipeline includes and pipeline roots in nf-lang"
```

---

### Task 2: `.pipeline-info` and pipeline checksum

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/pipeline/PipelineInfo.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/pipeline/PipelineInfoTest.groovy`

**Interfaces:**
- Produces: `class PipelineInfo`, with
  - fields `String name, version, workflowPath, source, repositoryUrl, revision, commitId, checksum, registryUrl`
  - `static PipelineInfo load(Path dir)`: null when the file is missing
  - `void save(Path dir)`
  - `static String checksum(Path dir)`
  - `boolean isModified(Path dir)`: `checksum(dir) != this.checksum`

- [ ] **Step 1: Write the failing tests:**
  - `save then load round-trips every field`
  - `load returns null when .pipeline-info is missing`
  - `checksum is stable and changes when a file changes`: write `main.nf`, take the checksum, edit `main.nf`, and assert the value changed.
  - `checksum ignores .pipeline-info and .git`
  - `ignores nested dirs with their own info file`: take the checksum, then create `modules/nf-core/fastqc/{main.nf,.module-info}` and `pipelines/x/y/{main.nf,.pipeline-info}`, and assert the checksum is unchanged. Adding `modules/nf-core/other/main.nf` *without* `.module-info` changes it.

- [ ] **Step 2:** Run `make test module=nextflow class=PipelineInfoTest`. Expected: FAIL.

- [ ] **Step 3: Implement.**
  - `load` and `save` use `java.util.Properties` on `RemoteModuleResolver.PIPELINE_INFO_FILE`. Null fields are not written.
  - `checksum` uses `Files.walkFileTree`:
    - `SKIP_SUBTREE` for `.git`.
    - `SKIP_SUBTREE` for any directory, other than `dir` itself, that contains `.module-info` or `.pipeline-info`.
    - Then use the same digest loop as `ModuleChecksum.computeDirectory`.

- [ ] **Step 4:** Run `make test module=nextflow class=PipelineInfoTest`. Expected: PASS.

- [ ] **Step 5:** Commit: `git commit -s -m "Add .pipeline-info and pipeline checksum"`.

---

### Task 3: Fetch a pinned commit with Git

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/pipeline/PipelineGitFetcher.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/pipeline/PipelineGitFetcherTest.groovy`

**Interfaces:**
- Produces: `class PipelineGitFetcher`, with `void fetch(String repositoryUrl, String revision, String commitId, Path target)`. Afterwards, `target` holds the exact tree at `commitId`, without `.git`. It throws `AbortOperationException` when the commit can't be found.

- [ ] **Step 1: Write the failing tests.** Use a JGit repo created in `@TempDir` with two commits and tag `1.0` on the first, and pass it as a `file://` URL.
  - `fetches the tag and removes .git`: `fetch(url, '1.0', <sha1>, target)`. Assert `target/main.nf` has the first commit's content and `!Files.exists(target.resolve('.git'))`.
  - `fetches a commit sha revision`: `fetch(url, <sha2>, <sha2>, target)`. Assert the second commit's content.
  - `fails when the revision no longer points to the pinned commit`: `fetch(url, '1.0', <sha2>, target)`. Expect `AbortOperationException` whose message contains `does not match the pinned commit`.

- [ ] **Step 2:** Run `make test module=nextflow class=PipelineGitFetcherTest`. Expected: FAIL.

- [ ] **Step 3: Implement with JGit.**
  - If `revision` is not a 40-hex SHA, clone with `setBranch(revision)` and `setDepth(1)`.
  - If it is a SHA, clone without a branch or depth, then check out `commitId`.
  - Then compare `HEAD` with `commitId`. On a mismatch, throw: `Revision '<revision>' of <url> does not match the pinned commit <commitId>`.
  - Finally delete `.git`.
  - Plain JGit is enough: the registry accepts only public repositories, so the `AssetManager` credentials and provider code is not needed. This is a known limit; switch to `AssetManager` if private registries arrive.

- [ ] **Step 4:** Run `make test module=nextflow class=PipelineGitFetcherTest`. Expected: PASS.

- [ ] **Step 5:** Commit: `git commit -s -m "Fetch pinned pipeline commits with JGit"`.

---

### Task 4: PipelineManager (resolve, install, list, remove)

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/pipeline/PipelineManager.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/pipeline/PipelineManagerTest.groovy`

**Interfaces:**
- Consumes: `PipelineInfo` (Task 2), `PipelineGitFetcher.fetch` (Task 3), and these `npr-client` methods: `RegistryClient.getPipeline(name): Pipeline` (with `getLatest()`, `getRepositoryUrl()`, `getWorkflowPath()`, `getRevoked()`) and `getPipelineRelease(name, version): PipelineRelease` (with `getVersion()`, `getGitRevision()`, `getCommitId()`, `getStatus()`).
- Produces: `class PipelineManager`, with
  - constructor `PipelineManager(Path baseDir, RegistryClient client, PipelineGitFetcher fetcher = new PipelineGitFetcher())`
  - `Path resolve(String name)`: auto-installs, and returns the main script path (`<dir>/<workflowPath>`)
  - `PipelineInfo install(String name, String version, boolean force)`: `version` null means latest
  - `List<PipelineInfo> list()`: only `<baseDir>/pipelines/*/*`
  - `boolean remove(String name, boolean force)`
  - `Path pipelineDir(String name)`

- [ ] **Step 1: Write the failing tests.** Use `Mock(RegistryClient)` and `Mock(PipelineGitFetcher)`; the fetcher mock writes `main.nf` into its `target` argument.
  - `installs latest when missing`: `getPipeline` returns latest `3.14.0`, revision `3.14.0`, commit `abc`. After `resolve('nf-core/rnaseq')`: the result is `pipelines/nf-core/rnaseq/main.nf`, `.pipeline-info` has `version=3.14.0`, `commitId=abc`, `source=git`, and a non-null `checksum`.
  - `uses installed pipeline without calling the registry`: pre-install, then `resolve` again. Assert `0 * client._`.
  - `warns and keeps a modified pipeline`: edit `main.nf`, call `resolve`. Assert the same path is returned, the file content is unchanged, and `0 * fetcher._`.
  - `install refuses a modified pipeline without force, replaces it with force`.
  - `install of the same version and checksum is a no-op` (`0 * fetcher._`).
  - `install -version replaces an unmodified different version`.
  - `failed install leaves no directory`: the fetcher throws. Assert `!Files.exists(pipelineDir)` and that a second `resolve` calls the fetcher again.
  - `fails on no published release / revoked`: if `latest` is null, expect the message to contain `has no published release`. If the pipeline is revoked or the release status is `REVOKED`, expect it to contain `is revoked`.
  - `rejects invalid pipeline names`: `'../x'`, `'nf-core/a/b'`, `'NF-core/x'`, `'x'`. Expect `AbortOperationException`, with no client calls.
  - `list returns installed pipelines`; `remove deletes, refuses modified without force`.

- [ ] **Step 2:** Run `make test module=nextflow class=PipelineManagerTest`. Expected: FAIL.

- [ ] **Step 3: Implement.**
  - **Name check:** `^[a-z0-9][a-z0-9._-]*/[a-z][a-z0-9._-]*$`.
  - **Install:**
    1. Fetch into `Files.createTempDirectory(baseDir/pipelines/<scope>, '.tmp-<name>-')`.
    2. Write `.pipeline-info` with its `checksum`.
    3. Delete the old directory if `force` or if it is unmodified.
    4. `Files.move(tmp, dir, ATOMIC_MOVE)`.
    5. On any error, delete `tmp`.
  - **Logging:** `log.warn1` for modified (`Pipeline <name> has local modifications (checksum mismatch)`). `log.info` for `Installed pipeline <name>@<version>`.

- [ ] **Step 4:** Run `make test module=nextflow class=PipelineManagerTest`. Expected: PASS.

- [ ] **Step 5:** Commit: `git commit -s -m "Add PipelineManager for pipeline install, list and remove"`.

---

### Task 5: Wire pipelines into compile and runtime

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/module/DefaultRemoteModuleResolver.groovy` (add `resolvePipeline`)
- Modify: `modules/nextflow/src/main/groovy/nextflow/script/IncludeDef.groovy:180-189` (`resolveRemoteModulePath`)
- Modify: `modules/nextflow/src/main/groovy/nextflow/module/ModuleResolver.groovy:75-79` (skip the "no registry origin" warning under a pipeline root)
- Test: `modules/nextflow/src/test/groovy/nextflow/script/PipelineCompositionTest.groovy`, `modules/nextflow/src/test/groovy/nextflow/module/ModuleResolverTest.groovy`

**Interfaces:**
- Consumes: `RemoteModuleResolver.isPipelineInclude`, `resolveBaseDir` (Task 1), `PipelineManager.resolve` (Task 4).

- [ ] **Step 1: Write the failing tests.**
  - `PipelineCompositionTest`, `should include a pipeline by registry name`:
    - Write `pipelines/test/pipe/{main.nf,.pipeline-info}`. The `.pipeline-info` has a checksum that matches the files, made with `PipelineInfo`.
    - Write a meta-pipeline `main.nf` with `include { workflow as PIPE } from 'test/pipe'` that calls `PIPE(...)`.
    - Run it as the existing composition tests do. Assert the output, and that no registry call happens (the default resolver with an unreachable `registry.url`).
  - `PipelineCompositionTest`, `should resolve a registry module inside a pipeline against the pipeline`:
    - `pipelines/test/pipe/main.nf` includes `{ HELLO } from 'test/hello'`.
    - Provide `pipelines/test/pipe/modules/test/hello/{main.nf,meta.yml}`, without `.module-info`.
    - Provide a *different* `modules/test/hello` at the project root that prints something else.
    - Assert that the pipeline's copy runs.
  - `ModuleResolverTest`, `does not warn for a module without .module-info inside a pipeline root`: capture the log and assert that `no registry origin` is absent.

- [ ] **Step 2:** Run `make test module=nextflow class=PipelineCompositionTest`. Expected: the new tests FAIL.

- [ ] **Step 3: Implement.**
  - `DefaultRemoteModuleResolver.resolvePipeline(name, baseDir)` builds a `PipelineManager(baseDir, RegistryClientFactory.forConfig(registryConfig))` with the same config lookup as `resolve`, and returns `manager.resolve(name)`. Errors are wrapped in `IllegalModulePathException`, as for modules.
  - `IncludeDef.resolveRemoteModulePath` calls `resolver.resolvePipeline` when `RemoteModuleResolver.isPipelineInclude(modules*.name)`.
  - `ModuleResolver`: skip the `NO_REMOTE_MODULE` warning when the module directory has a parent containing `.pipeline-info` below `baseDir`.

- [ ] **Step 4:** Run `make test module=nextflow`. Expected: all PASS.

- [ ] **Step 5:** Commit: `git commit -s -m "Auto-install registry pipelines on include"`.

---

### Task 6: `nextflow pipeline install|list|remove`

**Files:**
- Create: `modules/nf-cli-v1/src/main/groovy/nextflow/cli/CmdPipeline.groovy`
- Create: `modules/nf-cli-v1/src/main/groovy/nextflow/cli/pipeline/CmdPipelineInstall.groovy`, `CmdPipelineList.groovy`, `CmdPipelineRemove.groovy`
- Modify: `modules/nf-cli-v1/src/main/groovy/nextflow/cli/Launcher.groovy:117` (register `new CmdPipeline()`)
- Test: `modules/nf-cli-v1/src/test/groovy/nextflow/cli/pipeline/CmdPipelineInstallTest.groovy`, `CmdPipelineListTest.groovy`, `CmdPipelineRemoveTest.groovy`

**Interfaces:**
- Consumes: `PipelineManager` (Task 4).
- Produces:
  - `CmdPipeline`, with `NAME = 'pipeline'`, the same dispatch as `CmdModule`, and program name `nextflow pipeline`.
  - Each subcommand has `protected RegistryClient client` for tests, like `CmdModuleSearch`, and a `Path root` that defaults to the current directory.

- [ ] **Step 1: Write the failing tests.** Copy the setup of the matching `CmdModule*Test`.
  - Install: `install nf-core/rnaseq -version 3.14.0` calls `getPipelineRelease('nf-core/rnaseq','3.14.0')` and creates the directory. Without `-version` it calls `getPipeline`. On a modified copy without `-force` it throws `AbortOperationException` with a message that contains `-force`.
  - List: prints `nf-core/rnaseq  3.14.0  git  ok` (column text left to the implementer, but these four values must appear). `-output json` prints a JSON array with `name`, `version`, `source`, `status`.
  - Remove: deletes the directory. On a modified copy without `-force` it refuses.

- [ ] **Step 2:** Run `make test module=nf-cli-v1 class='nextflow.cli.pipeline.*'`. Expected: FAIL.

- [ ] **Step 3: Implement** the commands as thin wrappers over `PipelineManager`. Options:
  - `install`: `-version`, `-force`
  - `list`: `-output text|json`
  - `remove`: `-force`

- [ ] **Step 4:** Run `make test module=nf-cli-v1`. Expected: PASS.

- [ ] **Step 5:** Commit: `git commit -s -m "Add nextflow pipeline install, list and remove commands"`.

---

### Task 7: `nextflow pipeline search|view`

**Files:**
- Create: `modules/nf-cli-v1/src/main/groovy/nextflow/cli/pipeline/CmdPipelineSearch.groovy`, `CmdPipelineView.groovy`
- Modify: `CmdPipeline.groovy` (register them)
- Test: `CmdPipelineSearchTest.groovy`, `CmdPipelineViewTest.groovy`

**Interfaces:**
- Consumes: `RegistryClient.searchPipelines(query, limit): ListPipelinesResponse` (`getPipelines()` → `PipelineSummary` with `getName()`, `getLatestVersion()`, `getMetadata().getDescription()`), `getPipeline`, `getPipelineRelease` (with `getSchema(): Map`).

- [ ] **Step 1: Write the failing tests.**
  - Search: `search rnaseq -limit 5` calls `searchPipelines('rnaseq', 5)` and prints each name, version and description. `-output json` prints the list.
  - View: `view nf-core/rnaseq` prints `Name:`, `Version:`, `Repository:`, `Revision:`, `Commit:`, the registry URL `<registry>/pipelines/nf-core/rnaseq@<version>`, and the usage block below. When the schema has params, it also lists every param in `schema.$defs.*.properties` and `schema.properties`, with type and description.

```groovy
include { params as RnaseqParams ; workflow as NFCORE_RNASEQ } from 'nf-core/rnaseq'
```

    The alias is `NFCORE_` + the upper-cased name with `-` replaced by `_`, and the params type is the PascalCase name + `Params`. `-version` calls `getPipelineRelease`. `-output json` prints the release and pipeline fields.

- [ ] **Step 2:** Run `make test module=nf-cli-v1 class='nextflow.cli.pipeline.*'`. Expected: FAIL.

- [ ] **Step 3: Implement** following `CmdModuleSearch` and `CmdModuleView`. Use the same `-limit` default and `-output` validator as the module commands.

- [ ] **Step 4:** Run `make test module=nf-cli-v1`. Expected: PASS.

- [ ] **Step 5:** Commit: `git commit -s -m "Add nextflow pipeline search and view commands"`.

---

### Task 8: End-to-end check against the registry

- [ ] **Step 1:** `make compile`, then `./launch.sh pipeline search rnaseq` and `./launch.sh pipeline view nf-core/rnaseq`. Expected: results from `registry.nextflow.io`.
- [ ] **Step 2:** In a temp dir, run `./launch.sh pipeline install nf-core/rnaseq -version 3.26.0`. Expected: `pipelines/nf-core/rnaseq/.pipeline-info` with `commitId=e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`, then `./launch.sh pipeline list` shows `ok`.
- [ ] **Step 3:** Run `./launch.sh pipeline remove nf-core/rnaseq`. Expected: the directory is gone.
- [ ] **Step 4:** Once seqeralabs/nextflow-registry#439 is released, set `nprVersion` to the released version, remove the local `mavenLocal()` line, run `make test module=nextflow` and `make test module=nf-cli-v1`, and commit only `gradle.properties`: `git commit -s -m "Bump npr to <version> for pipeline registry methods"`.

## Not in this plan

- User docs for `nextflow pipeline`: there are no CLI docs for `nextflow module` in `docs/` yet either. Add them with the module docs.
- M1 (publish, bundle install): a separate plan, after seqeralabs/nextflow-registry#440.
