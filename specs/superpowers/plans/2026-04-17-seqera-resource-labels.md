# Seqera executor `process.resourceLabels` Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the `nf-seqera` executor honour `process.resourceLabels` by sending the config-level baseline as run labels and the per-task delta as Sched task labels.

**Architecture:** Cumulative Nextflow labels split into two scheduler scopes — config-level `process.resourceLabels` becomes `CreateRunRequest.labels`; the difference between `task.config.getResourceLabels()` and that baseline becomes `Task.labels`. The redundant `seqera.executor.labels` config option is removed.

**Tech Stack:** Groovy 4 / Java 21 toolchain, Gradle, Spock, `io.seqera:sched-client` (≥ 0.51.0 — must expose `Task.labels`), Nextflow extension-point plugin model.

**Spec:** `docs/superpowers/specs/2026-04-17-seqera-resource-labels-design.md`

**File map:**
- Modify `settings.gradle` — uncomment-style `includeBuild '../sched'` for dev
- Modify `plugins/nf-seqera/build.gradle` — bump `sched-client` version
- Modify `plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy` — add `withProcessResourceLabels`, add `delta`, remove `withUserLabels`
- Modify `plugins/nf-seqera/src/main/io/seqera/config/ExecutorOpts.groovy` — remove `labels` field
- Modify `plugins/nf-seqera/src/main/io/seqera/executor/SeqeraExecutor.groovy` — wire run labels + cache `runResourceLabels`
- Modify `plugins/nf-seqera/src/main/io/seqera/executor/SeqeraTaskHandler.groovy` — attach delta to `Task.labels`
- Modify `plugins/nf-seqera/src/test/io/seqera/executor/LabelsTest.groovy`
- Modify `plugins/nf-seqera/src/test/io/seqera/config/ExecutorOptsTest.groovy`
- Modify `plugins/nf-seqera/src/test/io/seqera/executor/SeqeraTaskHandlerTest.groovy`
- Modify `docs/reference/process.md` — add Seqera executor to support list
- Modify `plugins/nf-seqera/changelog.txt` — entry
- Modify `plugins/nf-seqera/VERSION` — bump to `0.18.0`

---

### Task 1: Bump `sched-client` and wire `includeBuild '../sched'`

**Files:**
- Modify: `plugins/nf-seqera/build.gradle:54`
- Modify: `settings.gradle:17-19`

The local `~/Projects/sched` checkout is at `0.51.0`; that is the version exposing `Task.labels`. Until 0.51.0 is published to the Seqera Maven repo, we use a Gradle composite build to substitute the dependency from the local checkout.

- [ ] **Step 1: Bump sched-client version**

Edit `plugins/nf-seqera/build.gradle:54`:

```gradle
    api 'io.seqera:sched-client:0.51.0'
```

- [ ] **Step 2: Add `includeBuild '../sched'` block to `settings.gradle`**

Replace the commented `pluginManagement` block at `settings.gradle:17-19` with both blocks (keep the existing comment, add a new one for sched as dev-only opt-in):

```gradle
// pluginManagement {
//     includeBuild '../nextflow-plugin-gradle'
// }

// For local development against an unpublished sched-client, uncomment:
// includeBuild '../sched'
includeBuild '../sched'
```

(The uncommented `includeBuild '../sched'` line is required for the build to resolve `sched-client:0.51.0` until the artifact is published. The commented hint stays for future reference.)

- [ ] **Step 3: Verify the build resolves**

Run: `./gradlew :plugins:nf-seqera:compileGroovy`
Expected: BUILD SUCCESSFUL. If it fails with a missing `sched-client:0.51.0`, confirm `~/Projects/sched/VERSION` contains `0.51.0` and that `~/Projects/sched/sched-client` builds locally (`cd ~/Projects/sched && ./gradlew :sched-client:assemble`).

- [ ] **Step 4: Commit**

```bash
git add settings.gradle plugins/nf-seqera/build.gradle
git commit -s -m "build(nf-seqera): bump sched-client to 0.51.0 via includeBuild"
```

---

### Task 2: Add `Labels.withProcessResourceLabels` (TDD)

**Files:**
- Modify: `plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy`
- Modify: `plugins/nf-seqera/src/test/io/seqera/executor/LabelsTest.groovy`

- [ ] **Step 1: Write the failing tests**

Append to `plugins/nf-seqera/src/test/io/seqera/executor/LabelsTest.groovy` (before the closing `}`):

```groovy
    def 'should add process resource labels coercing values to string'() {
        when:
        def labels = new Labels()
                .withProcessResourceLabels([team: 'genomics', priority: 7, retain: true])

        then:
        labels.entries['team'] == 'genomics'
        labels.entries['priority'] == '7'
        labels.entries['retain'] == 'true'
    }

    def 'should ignore null or empty process resource labels'() {
        when:
        def a = new Labels().withProcessResourceLabels(null)
        def b = new Labels().withProcessResourceLabels([:])

        then:
        a.entries.isEmpty()
        b.entries.isEmpty()
    }

    def 'should let process resource labels override workflow metadata on key collision'() {
        given:
        def workflow = Mock(WorkflowMetadata) {
            getProjectName() >> 'hello'
            getRunName() >> 'happy_turing'
            getSessionId() >> UUID.randomUUID()
            isResume() >> false
            getManifest() >> new Manifest([:])
        }

        when:
        def labels = new Labels()
                .withWorkflowMetadata(workflow)
                .withProcessResourceLabels(['nextflow.io/runName': 'custom', team: 'a'])

        then:
        labels.entries['nextflow.io/runName'] == 'custom'
        labels.entries['team'] == 'a'
        labels.entries['nextflow.io/projectName'] == 'hello'
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `./gradlew :plugins:nf-seqera:test --tests 'io.seqera.executor.LabelsTest' -i`
Expected: FAIL — `MissingMethodException: No signature of method ... withProcessResourceLabels`.

- [ ] **Step 3: Implement `withProcessResourceLabels`**

Edit `plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy`. After the `withUserLabels` method (which will be removed in Task 4), add:

```groovy
    /**
     * Add config-level {@code process.resourceLabels}. Values are coerced to
     * string via {@link String#valueOf} to satisfy the scheduler API typing.
     */
    Labels withProcessResourceLabels(Map<String,?> map) {
        if( !map ) return this
        map.each { k, v -> entries.put(k.toString(), String.valueOf(v)) }
        return this
    }
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `./gradlew :plugins:nf-seqera:test --tests 'io.seqera.executor.LabelsTest' -i`
Expected: PASS for the three new tests; existing tests still PASS.

- [ ] **Step 5: Commit**

```bash
git add plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy \
        plugins/nf-seqera/src/test/io/seqera/executor/LabelsTest.groovy
git commit -s -m "feat(nf-seqera): add Labels.withProcessResourceLabels"
```

---

### Task 3: Add `Labels.delta` and `Labels.toStringMap` helpers (TDD)

**Files:**
- Modify: `plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy`
- Modify: `plugins/nf-seqera/src/test/io/seqera/executor/LabelsTest.groovy`

These two static helpers compute the per-task delta and coerce arbitrary `Map<String,?>` values to strings — used by both the executor (to cache the run baseline) and the task handler (to compute the delta).

- [ ] **Step 1: Write the failing tests**

Append to `plugins/nf-seqera/src/test/io/seqera/executor/LabelsTest.groovy` (before the closing `}`):

```groovy
    def 'should coerce map values to strings'() {
        expect:
        Labels.toStringMap(null) == [:]
        Labels.toStringMap([:]) == [:]
        Labels.toStringMap([a: 1, b: 'x', c: true]) == [a: '1', b: 'x', c: 'true']
    }

    def 'should compute null delta when task labels are empty'() {
        expect:
        Labels.delta(null, [team: 'a']) == null
        Labels.delta([:], [team: 'a']) == null
    }

    def 'should return full task labels when run labels are empty'() {
        expect:
        Labels.delta([team: 'a', region: 'us'], null) == [team: 'a', region: 'us']
        Labels.delta([team: 'a', region: 'us'], [:]) == [team: 'a', region: 'us']
    }

    def 'should keep only differing or missing keys in delta'() {
        expect:
        Labels.delta([team: 'a', region: 'us'], [team: 'a']) == [region: 'us']
        Labels.delta([team: 'b'], [team: 'a']) == [team: 'b']
        Labels.delta([team: 'a', region: 'us'], [team: 'a', region: 'us']) == null
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `./gradlew :plugins:nf-seqera:test --tests 'io.seqera.executor.LabelsTest' -i`
Expected: FAIL — `MissingMethodException: ... toStringMap` and `... delta`.

- [ ] **Step 3: Implement the helpers**

Edit `plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy`. Inside the class, after the existing `runId` method, add:

```groovy
    /**
     * Coerce arbitrary map values to strings via {@link String#valueOf}.
     * Returns an empty map for null/empty input.
     */
    static Map<String,String> toStringMap(Map<String,?> map) {
        if( !map ) return Collections.<String,String>emptyMap()
        final result = new LinkedHashMap<String,String>(map.size())
        map.each { k, v -> result.put(k.toString(), String.valueOf(v)) }
        return result
    }

    /**
     * Return the entries of {@code task} that are missing from {@code run}
     * or have a different value. Returns {@code null} if the resulting
     * map would be empty (so callers can omit the field).
     */
    static Map<String,String> delta(Map<String,String> task, Map<String,String> run) {
        if( !task ) return null
        final result = new LinkedHashMap<String,String>()
        task.each { k, v ->
            if( run == null || !run.containsKey(k) || run.get(k) != v )
                result.put(k, v)
        }
        return result.isEmpty() ? null : result
    }
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `./gradlew :plugins:nf-seqera:test --tests 'io.seqera.executor.LabelsTest' -i`
Expected: PASS for all new tests; existing tests still PASS.

- [ ] **Step 5: Commit**

```bash
git add plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy \
        plugins/nf-seqera/src/test/io/seqera/executor/LabelsTest.groovy
git commit -s -m "feat(nf-seqera): add Labels.toStringMap and Labels.delta helpers"
```

---

### Task 4: Remove `seqera.executor.labels` config option

**Files:**
- Modify: `plugins/nf-seqera/src/main/io/seqera/config/ExecutorOpts.groovy:74-79,130,168-170`
- Modify: `plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy:81-88`
- Modify: `plugins/nf-seqera/src/main/io/seqera/executor/SeqeraExecutor.groovy:120`
- Modify: `plugins/nf-seqera/src/test/io/seqera/config/ExecutorOptsTest.groovy:131-165`
- Modify: `plugins/nf-seqera/src/test/io/seqera/executor/LabelsTest.groovy:127-150,192-199`

The user-facing `seqera.executor.labels` option is replaced by the standard Nextflow `process.resourceLabels` directive.

- [ ] **Step 1: Remove the field and getter from `ExecutorOpts`**

Edit `plugins/nf-seqera/src/main/io/seqera/config/ExecutorOpts.groovy`:

Remove lines 74-79 (the `@ConfigOption` block and `final Map<String, String> labels` field):

```groovy
    @ConfigOption
    @Description("""
        Custom labels to apply to AWS resources for cost tracking and resource organization.
        Labels are propagated to ECS tasks, capacity providers, and EC2 instances.
    """)
    final Map<String, String> labels
```

Remove the assignment in the constructor (around line 129-130):

```groovy
        // labels for cost tracking
        this.labels = opts.labels as Map<String, String>
```

Remove the getter (around line 168-170):

```groovy
    Map<String, String> getLabels() {
        return labels
    }
```

- [ ] **Step 2: Remove the `withUserLabels` method from `Labels`**

Edit `plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy`. Delete the entire `withUserLabels` method (lines 81-88):

```groovy
    /**
     * Add user-configured labels. These take precedence over implicit labels.
     */
    Labels withUserLabels(Map<String,String> labels) {
        if( labels )
            entries.putAll(labels)
        return this
    }
```

- [ ] **Step 3: Remove the `withUserLabels` call site in `SeqeraExecutor.createRun()`**

Edit `plugins/nf-seqera/src/main/io/seqera/executor/SeqeraExecutor.groovy`. Delete the line:

```groovy
        labels.withUserLabels(seqeraConfig.labels)
```

- [ ] **Step 4: Remove obsolete tests**

Edit `plugins/nf-seqera/src/test/io/seqera/config/ExecutorOptsTest.groovy`. Delete the three tests (lines 131-165): `'should create config with labels'`, `'should handle null labels'`, `'should handle empty labels'`.

Edit `plugins/nf-seqera/src/test/io/seqera/executor/LabelsTest.groovy`. Delete the two tests: `'should allow user labels to override implicit labels'` (lines 127-150) and `'should handle null user labels'` (lines 192-199).

- [ ] **Step 5: Compile and run tests**

Run: `./gradlew :plugins:nf-seqera:compileGroovy :plugins:nf-seqera:test`
Expected: BUILD SUCCESSFUL; all remaining tests PASS. If the compiler complains about a stray reference to `seqeraConfig.labels` or `withUserLabels`, grep for and remove them: `rg "seqeraConfig\.labels|withUserLabels" plugins/nf-seqera`.

- [ ] **Step 6: Commit**

```bash
git add plugins/nf-seqera/src/main/io/seqera/config/ExecutorOpts.groovy \
        plugins/nf-seqera/src/main/io/seqera/executor/Labels.groovy \
        plugins/nf-seqera/src/main/io/seqera/executor/SeqeraExecutor.groovy \
        plugins/nf-seqera/src/test/io/seqera/config/ExecutorOptsTest.groovy \
        plugins/nf-seqera/src/test/io/seqera/executor/LabelsTest.groovy
git commit -s -m "refactor(nf-seqera)!: remove seqera.executor.labels in favour of process.resourceLabels"
```

---

### Task 5: Wire `process.resourceLabels` into `SeqeraExecutor.createRun()` and expose `runResourceLabels` (TDD)

**Files:**
- Modify: `plugins/nf-seqera/src/main/io/seqera/executor/SeqeraExecutor.groovy`
- Modify: `plugins/nf-seqera/src/test/io/seqera/executor/SeqeraExecutorTest.groovy`

The executor reads the config-level `process.resourceLabels` map once at run creation, attaches it to the run labels via `Labels.withProcessResourceLabels`, and caches the coerced map so task handlers can compute deltas.

- [ ] **Step 1: Write the failing test**

Append to `plugins/nf-seqera/src/test/io/seqera/executor/SeqeraExecutorTest.groovy` (before the final closing `}`):

```groovy
    def 'should expose run resource labels coerced from config-level process.resourceLabels'() {
        given:
        def executor = new SeqeraExecutor()
        executor.@session = Mock(Session) {
            getConfig() >> [process: [resourceLabels: [team: 'a', priority: 7]]]
        }

        when:
        executor.computeRunResourceLabels()

        then:
        executor.runResourceLabels == [team: 'a', priority: '7']
    }

    def 'should yield empty run resource labels when process.resourceLabels is absent'() {
        given:
        def executor = new SeqeraExecutor()
        executor.@session = Mock(Session) {
            getConfig() >> [:]
        }

        when:
        executor.computeRunResourceLabels()

        then:
        executor.runResourceLabels == [:]
    }
```

(`Session` is already imported at line 21; if not, add the import.)

- [ ] **Step 2: Run the test to verify it fails**

Run: `./gradlew :plugins:nf-seqera:test --tests 'io.seqera.executor.SeqeraExecutorTest' -i`
Expected: FAIL — `computeRunResourceLabels` / `runResourceLabels` don't exist.

- [ ] **Step 3: Implement on `SeqeraExecutor`**

Edit `plugins/nf-seqera/src/main/io/seqera/executor/SeqeraExecutor.groovy`.

Add a private field near the other private fields (after `runId` at line 65):

```groovy
    private volatile Map<String,String> runResourceLabels = Collections.<String,String>emptyMap()
```

Add a method to compute the run resource labels (place near other protected/package methods, e.g. before `createRun()` at line 110):

```groovy
    @groovy.transform.PackageScope
    void computeRunResourceLabels() {
        final processMap = session.config.process as Map
        final raw = processMap?.get('resourceLabels') as Map<String,?>
        this.runResourceLabels = Labels.toStringMap(raw)
    }
```

Add the public getter (after `getRunId()` around line 204):

```groovy
    Map<String,String> getRunResourceLabels() {
        return runResourceLabels
    }
```

Wire it into `createRun()`. Replace the labels-building block at `SeqeraExecutor.groovy:117-120` (after the deletion in Task 4 it should look like the first three lines below) with:

```groovy
        computeRunResourceLabels()
        final labels = new Labels()
        if( seqeraConfig.autoLabels )
            labels.withWorkflowMetadata(session.workflowMetadata)
        labels.withProcessResourceLabels(runResourceLabels)
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `./gradlew :plugins:nf-seqera:test --tests 'io.seqera.executor.SeqeraExecutorTest' -i`
Expected: PASS for both new tests; existing tests still PASS.

- [ ] **Step 5: Commit**

```bash
git add plugins/nf-seqera/src/main/io/seqera/executor/SeqeraExecutor.groovy \
        plugins/nf-seqera/src/test/io/seqera/executor/SeqeraExecutorTest.groovy
git commit -s -m "feat(nf-seqera): attach process.resourceLabels to Sched run labels"
```

---

### Task 6: Send per-task delta on `Task.labels` from `SeqeraTaskHandler.submit()` (TDD)

**Files:**
- Modify: `plugins/nf-seqera/src/main/io/seqera/executor/SeqeraTaskHandler.groovy`
- Modify: `plugins/nf-seqera/src/test/io/seqera/executor/SeqeraTaskHandlerTest.groovy`

Capture the `Task` passed to the batch submitter and assert its `labels` field reflects the delta between the task's `getResourceLabels()` and the executor's `runResourceLabels`.

- [ ] **Step 1: Write the failing tests**

Append to `plugins/nf-seqera/src/test/io/seqera/executor/SeqeraTaskHandlerTest.groovy` (before the final closing `}`):

```groovy
    def 'submit attaches Task.labels containing only the per-task delta'() {
        given:
        Task captured = null
        def batchSubmitter = Mock(SeqeraBatchSubmitter) {
            submit(_, _) >> { args -> captured = args[1] as Task }
        }
        def taskConfig = Mock(TaskConfig) {
            getCpus() >> 2
            getMemory() >> MemoryUnit.of('1 GB')
            getAccelerator() >> null
            getResourceLabels() >> [team: 'a', region: 'us-east-1']
            getResourceLimit('memory') >> null
            getResourceLimit('cpus') >> null
            getDisk() >> null
        }
        def taskRun = Mock(TaskRun) {
            getConfig() >> taskConfig
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getWorkDirStr() >> '/work/ab/cd1234'
            getContainer() >> 'docker.io/library/alpine:3'
            getContainerPlatform() >> 'linux/amd64'
            getId() >> TaskId.of(1)
            getHash() >> HashCode.fromInt(1)
            lazyName() >> 'sample_task'
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
            getBatchSubmitter() >> batchSubmitter
            getSeqeraConfig() >> Mock(ExecutorOpts) {
                getMachineRequirement() >> Mock(io.seqera.config.MachineRequirementOpts)
                getTaskEnvironment() >> [:]
            }
            getRunResourceLabels() >> [team: 'a']
            ensureRunCreated() >> {}
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            fusionEnabled() >> true
            fusionLauncher() >> Mock(nextflow.fusion.FusionScriptLauncher) {
                fusionEnv() >> [:]
            }
            fusionSubmitCli() >> ['/bin/sh', '-c', 'true']
            fusionConfig() >> Mock(nextflow.fusion.FusionConfig) {
                snapshotsEnabled() >> false
            }
        }

        when:
        handler.submit()

        then:
        captured != null
        captured.getLabels() == [region: 'us-east-1']
    }

    def 'submit leaves Task.labels unset when the task labels equal the run baseline'() {
        given:
        Task captured = null
        def batchSubmitter = Mock(SeqeraBatchSubmitter) {
            submit(_, _) >> { args -> captured = args[1] as Task }
        }
        def taskConfig = Mock(TaskConfig) {
            getCpus() >> 2
            getMemory() >> MemoryUnit.of('1 GB')
            getAccelerator() >> null
            getResourceLabels() >> [team: 'a']
            getResourceLimit('memory') >> null
            getResourceLimit('cpus') >> null
            getDisk() >> null
        }
        def taskRun = Mock(TaskRun) {
            getConfig() >> taskConfig
            getWorkDir() >> Paths.get('/work/ab/cd1234')
            getWorkDirStr() >> '/work/ab/cd1234'
            getContainer() >> 'docker.io/library/alpine:3'
            getContainerPlatform() >> 'linux/amd64'
            getId() >> TaskId.of(1)
            getHash() >> HashCode.fromInt(1)
            lazyName() >> 'sample_task'
        }
        def executor = Mock(SeqeraExecutor) {
            getClient() >> Mock(SchedClient)
            getBatchSubmitter() >> batchSubmitter
            getSeqeraConfig() >> Mock(ExecutorOpts) {
                getMachineRequirement() >> Mock(io.seqera.config.MachineRequirementOpts)
                getTaskEnvironment() >> [:]
            }
            getRunResourceLabels() >> [team: 'a']
            ensureRunCreated() >> {}
        }
        def handler = Spy(new SeqeraTaskHandler(taskRun, executor)) {
            fusionEnabled() >> true
            fusionLauncher() >> Mock(nextflow.fusion.FusionScriptLauncher) {
                fusionEnv() >> [:]
            }
            fusionSubmitCli() >> ['/bin/sh', '-c', 'true']
            fusionConfig() >> Mock(nextflow.fusion.FusionConfig) {
                snapshotsEnabled() >> false
            }
        }

        when:
        handler.submit()

        then:
        captured != null
        captured.getLabels() == null
    }
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `./gradlew :plugins:nf-seqera:test --tests 'io.seqera.executor.SeqeraTaskHandlerTest' -i`
Expected: FAIL — assertions on `captured.getLabels()` fail because submit() does not set them.

- [ ] **Step 3: Wire the delta into `submit()`**

Edit `plugins/nf-seqera/src/main/io/seqera/executor/SeqeraTaskHandler.groovy`. After the `final schedTask = new Task() ... .nextflow(...)` block ending around line 140, before the `log.debug` call at line 141, insert:

```groovy
        // attach per-task resource labels delta (over run-level baseline)
        final taskLabels = Labels.toStringMap(task.config.getResourceLabels())
        final delta = Labels.delta(taskLabels, executor.runResourceLabels)
        if( delta )
            schedTask.labels(delta)
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `./gradlew :plugins:nf-seqera:test --tests 'io.seqera.executor.SeqeraTaskHandlerTest' -i`
Expected: PASS for both new tests; existing tests still PASS.

- [ ] **Step 5: Run the full plugin test suite**

Run: `./gradlew :plugins:nf-seqera:test`
Expected: BUILD SUCCESSFUL; no regressions.

- [ ] **Step 6: Commit**

```bash
git add plugins/nf-seqera/src/main/io/seqera/executor/SeqeraTaskHandler.groovy \
        plugins/nf-seqera/src/test/io/seqera/executor/SeqeraTaskHandlerTest.groovy
git commit -s -m "feat(nf-seqera): send per-task resourceLabels delta on Sched task"
```

---

### Task 7: Docs, changelog, and version bump

**Files:**
- Modify: `docs/reference/process.md:1388-1393`
- Modify: `plugins/nf-seqera/changelog.txt`
- Modify: `plugins/nf-seqera/VERSION`

- [ ] **Step 1: Update docs**

Edit `docs/reference/process.md`. Replace the executor support list at lines 1388-1393:

```markdown
Resource labels are currently supported by the following executors:

- {ref}`awsbatch-executor`
- {ref}`azurebatch-executor`
- {ref}`google-batch-executor`
- {ref}`k8s-executor`
- {ref}`seqera-executor`
```

(If `seqera-executor` is not a defined ref, drop the `{ref}` wrapper and write `Seqera executor` as plain text.)

- [ ] **Step 2: Update plugin changelog**

Edit `plugins/nf-seqera/changelog.txt`. Add a new entry at the top, above the `0.17.0` block:

```
0.18.0 - <today's date>
- Support process.resourceLabels: config-level labels attached to Sched run, per-task delta attached to Sched task
- Remove seqera.executor.labels config option (use process.resourceLabels instead)
- Bump sched-client@0.51.0
```

- [ ] **Step 3: Bump plugin VERSION**

Edit `plugins/nf-seqera/VERSION`:

```
0.18.0
```

- [ ] **Step 4: Verify everything builds**

Run: `./gradlew :plugins:nf-seqera:check`
Expected: BUILD SUCCESSFUL.

- [ ] **Step 5: Commit**

```bash
git add docs/reference/process.md plugins/nf-seqera/changelog.txt plugins/nf-seqera/VERSION
git commit -s -m "docs(nf-seqera): document resourceLabels support and bump to 0.18.0"
```

---

## Self-review checklist (executed)

- **Spec coverage:** every section of `2026-04-17-seqera-resource-labels-design.md` maps to a task — sched-client bump (Task 1), `withProcessResourceLabels` (Task 2), `delta` + `toStringMap` (Task 3), removal of `seqera.executor.labels` (Task 4), run-level wiring + `runResourceLabels` (Task 5), per-task delta on `Task.labels` (Task 6), docs / changelog / VERSION (Task 7).
- **Placeholder scan:** no TBDs, no "implement later", every code step has the actual code.
- **Type consistency:** `Labels.toStringMap(Map<String,?>)` and `Labels.delta(Map<String,String>, Map<String,String>)` referenced consistently in Tasks 3, 5, 6; `runResourceLabels` field, `computeRunResourceLabels()` method, and `getRunResourceLabels()` getter consistent across Tasks 5 and 6.
