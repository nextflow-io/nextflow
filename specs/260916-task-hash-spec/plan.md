# Spec-driven versioned task hasher — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the family of `TaskHasher` subclasses with one interpreter (`BaseTaskHasher`) driven by data (`TaskHashSpec`), so hash versions and the global-cache hasher become spec values, and per-key digests become available for cross-version cache-miss diffing.

**Architecture:** A `TaskHashSpec` holds an ordered list of `KeyBinding` (a `HashKey` from a closed enum, plus a `Contributor` that emits 0..N values), an `EncodingRules` pair configuring `HashBuilder`, and a fingerprint derived from the canonical form of all of it. `BaseTaskHasher` folds every emitted value in order to produce the hash, or digests each binding separately to explain it. Master's existing `TaskHasher` is left untouched and becomes the byte-exactness test oracle.

**Tech Stack:** Groovy 4 / Java 17, Gradle, Spock, Guava `Hashing`/`Hasher`.

**Spec:** `specs/260916-task-hash-spec/spec.md`

## Status — executed, then restructured (2026-09-23)

**All 11 tasks were executed** via subagent-driven development on branch
`task-hash-spec-prototype`, each with a task review and, where findings arose, fix rounds. The branch
was then **rebased onto master** and restructured. Read the tasks below as the record of how it was
built, not as work still to do.

What the rebase and restructure changed:

| Plan task | Fate |
| --- | --- |
| 1–7 | Applied unchanged. |
| 8 (SPI + `TaskProcessor` wiring) | **Dropped.** Master landed `TaskHasherFactory` and `createTaskHasher()` independently, including the same `volatile` factory-list fix this task had made. Two commits were skipped during the rebase. |
| 9 (named `-dump-hashes`) | Hasher-side kept; its `TaskProcessor` edits dropped with task 8. |
| 10 (validation pipeline) | Applied; files live untracked under `specs/260916-task-hash-spec/pipeline/`. |
| 11 (`global/v1`) | Never dispatched — `GlobalTaskHasher` and `FileIdentityStrategy` exist only on seqeralabs#47, so there is nothing on this branch to modify. |

Two commits were then added on top, superseding parts of the design below:

- **Opt-in through master's hook.** `BaseTaskHasher` → `SpecTaskHasher extends TaskHasher`;
  `SpecTaskHasherFactory` abstains unless `NXF_TASK_HASH_VER` names a version; registered in core's
  `META-INF/extensions.idx`. `TaskProcessor` is byte-identical to master. Task 8's `TaskHashSpecFactory`
  is deleted.
- **Specs as JSON resources** under `nextflow/processor/hash/`, resolved through `ContributorRegistry`
  and validated by `TaskHashSpecLoader`. Task 5's and task 6's hard-coded `StdSpecs` binding lists are
  replaced by loaded files; the specs themselves are unchanged, so the oracle still passes.

Consequently the **Files** and **Interfaces** blocks in tasks 5, 6, 8 and 9 below are stale in their
class names and their `TaskProcessor` claims. The *reasoning* in each task — especially the four
byte-exactness failure modes in task 4 and the oracle design in task 5 — is unchanged and is why the
result can be trusted.

See `spec.md` for the current design, and `runbook.md` for validation results against genuine
releases.

## Global Constraints

- **Byte-exact reproduction is the acceptance test.** Any spec representing an existing hasher must produce an identical hash for the same task. No cleanup — no normalising absent keys, no reordering, no unifying encodings — inside a spec that reproduces existing behaviour.
- **`explain()` is strictly additive.** Nothing about the final hash bytes may depend on the per-key digest feature existing.
- **Default behaviour is unchanged.** With no env var and no plugin, a run uses `std/v4` and must hash identically to master.
- **All new files carry the Apache 2.0 licence header** used throughout the repo (`Copyright 2013-2026, Seqera Labs`).
- **New Groovy code is `@CompileStatic`** and follows the repo style: annotations on their own line, no single-line `if` bodies.
- **Commits are signed off** (`git commit -s`).
- `HashMode` is resolved at compute time from `task.processor.getConfig().getHashMode()`. It is never a spec field.

---

### Task 1: HashBuilder encoding flags and EncodingRules

The two `HashBuilder` behaviours that #6679 changed become switchable, defaulting to current behaviour so master is unaffected.

**Files:**
- Modify: `modules/nf-commons/src/main/nextflow/util/HashBuilder.java`
- Create: `modules/nf-commons/src/main/nextflow/util/EncodingRules.java`
- Test: `modules/nf-commons/src/test/nextflow/util/EncodingRulesTest.groovy`

**Interfaces:**
- Produces: `HashBuilder.withOrderIndependentMaps(boolean)`, `HashBuilder.withCacheFunnelFirst(boolean)`, both returning `HashBuilder`. `EncodingRules.RECORD_TYPES`, `EncodingRules.LEGACY`, `EncodingRules.apply(HashBuilder)`, `EncodingRules.canonicalForm()`.

- [ ] **Step 1: Write the failing test**

```groovy
// modules/nf-commons/src/test/nextflow/util/EncodingRulesTest.groovy
package nextflow.util

import com.google.common.hash.Hashing
import spock.lang.Specification

class EncodingRulesTest extends Specification {

    def 'legacy rules hash a map by values, record-types rules by entries'() {
        given:
        def map = [alpha: 'one', beta: 'two']

        when:
        def legacy = EncodingRules.LEGACY
            .apply(new HashBuilder().withHasher(Hashing.murmur3_128().newHasher()))
            .with(map)
            .build()
        and:
        def recordTypes = EncodingRules.RECORD_TYPES
            .apply(new HashBuilder().withHasher(Hashing.murmur3_128().newHasher()))
            .with(map)
            .build()

        then:
        legacy != recordTypes
    }

    def 'legacy rules hash only the values, so renaming a key does not change the hash'() {
        given:
        def one = [alpha: 'x']
        def two = [renamed: 'x']

        expect:
        EncodingRules.LEGACY.apply(new HashBuilder()).with(one).build() ==
            EncodingRules.LEGACY.apply(new HashBuilder()).with(two).build()
        and:
        EncodingRules.RECORD_TYPES.apply(new HashBuilder()).with(one).build() !=
            EncodingRules.RECORD_TYPES.apply(new HashBuilder()).with(two).build()
    }

    def 'default HashBuilder behaviour is unchanged'() {
        given:
        def map = [alpha: 'one', beta: 'two']

        expect:
        new HashBuilder().with(map).build() ==
            EncodingRules.RECORD_TYPES.apply(new HashBuilder()).with(map).build()
    }

    def 'canonical form is stable and distinguishes the rule sets'() {
        expect:
        EncodingRules.LEGACY.canonicalForm() == 'orderIndependentMaps=false;cacheFunnelFirst=false'
        and:
        EncodingRules.RECORD_TYPES.canonicalForm() == 'orderIndependentMaps=true;cacheFunnelFirst=true'
    }
}
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `./gradlew :nf-commons:test --tests "nextflow.util.EncodingRulesTest"`
Expected: FAIL — `unable to resolve class EncodingRules`.

- [ ] **Step 3: Add the flags to HashBuilder**

Add two fields beside the existing `mode` and `basePath` fields (around line 84), defaulting to current behaviour:

```java
    private boolean orderIndependentMaps = true;

    private boolean cacheFunnelFirst = true;
```

Add two builder methods beside `withBasePath`:

```java
    public HashBuilder withOrderIndependentMaps(boolean value) {
        this.orderIndependentMaps = value;
        return this;
    }

    public HashBuilder withCacheFunnelFirst(boolean value) {
        this.cacheFunnelFirst = value;
        return this;
    }
```

In `with(Object value)`, replace the `CacheFunnel` branch (currently immediately after the `Object[]` branch) with a guarded one, and make the `Map` branch conditional:

```java
        else if( cacheFunnelFirst && value instanceof CacheFunnel )
            ((CacheFunnel)value).funnel(hasher, mode);

        else if( value instanceof Map ) {
            if( orderIndependentMaps ) {
                hashUnorderedCollection(hasher, ((Map) value).entrySet(), mode);
            }
            else {
                // note: pre-#6679 behaviour — the map contributes its values only
                for( Object item : ((Map)value).values() ) {
                    with(item);
                }
            }
        }
```

Then add the late `CacheFunnel` branch immediately *after* the existing `SerializableMarker` branch and *before* the `Enum` branch, which is where it sat before #6679:

```java
        else if( !cacheFunnelFirst && value instanceof CacheFunnel )
            ((CacheFunnel)value).funnel(hasher, mode);
```

- [ ] **Step 4: Create EncodingRules**

```java
// modules/nf-commons/src/main/nextflow/util/EncodingRules.java
package nextflow.util;

/**
 * The per-key encoding half of a task hash version: how {@link HashBuilder} turns
 * any object into bytes. Frozen per {@code TaskHashSpec} — never mutated in place,
 * because historical rule sets must keep producing historical hashes.
 */
public class EncodingRules {

    /** Behaviour introduced by #6679 "Record types" (2026-03-09) — current master. */
    public static final EncodingRules RECORD_TYPES = new EncodingRules(true, true);

    /** Behaviour before #6679. */
    public static final EncodingRules LEGACY = new EncodingRules(false, false);

    private final boolean orderIndependentMaps;

    private final boolean cacheFunnelFirst;

    public EncodingRules(boolean orderIndependentMaps, boolean cacheFunnelFirst) {
        this.orderIndependentMaps = orderIndependentMaps;
        this.cacheFunnelFirst = cacheFunnelFirst;
    }

    public HashBuilder apply(HashBuilder builder) {
        return builder
            .withOrderIndependentMaps(orderIndependentMaps)
            .withCacheFunnelFirst(cacheFunnelFirst);
    }

    /** Stable text form, contributing to the spec fingerprint. Never reformat this. */
    public String canonicalForm() {
        return "orderIndependentMaps=" + orderIndependentMaps
            + ";cacheFunnelFirst=" + cacheFunnelFirst;
    }
}
```

- [ ] **Step 5: Run the tests and confirm they pass**

Run: `./gradlew :nf-commons:test --tests "nextflow.util.EncodingRulesTest" --tests "nextflow.util.HashBuilderTest"`
Expected: PASS, including the pre-existing `HashBuilderTest` — this is the proof the defaults are unchanged.

- [ ] **Step 6: Commit**

```bash
git add modules/nf-commons/src/main/nextflow/util/HashBuilder.java \
        modules/nf-commons/src/main/nextflow/util/EncodingRules.java \
        modules/nf-commons/src/test/nextflow/util/EncodingRulesTest.groovy
git commit -s -m "Add switchable HashBuilder encoding rules"
```

---

### Task 2: Key vocabulary and contributor types

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/hash/HashKey.groovy`
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/hash/Contributor.groovy`
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/hash/KeyBinding.groovy`
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/hash/HashContext.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/processor/hash/HashContextTest.groovy`

**Interfaces:**
- Consumes: nothing.
- Produces: `enum HashKey` (13 constants); `interface Contributor { String canonicalName(); List<Object> emit(HashContext ctx) }`; `class KeyBinding { HashKey key; Contributor contributor }`; `class HashContext` with `getTask()`, `getProcessor()`, `getSession()`, `globalVars()` returning `Map<String,Object>`, `binEntries()` returning `List<Path>`.

`HashContext` deliberately borrows `getTaskGlobalVars()` and `getTaskBinEntries(String)` from the existing `TaskHasher` rather than duplicating them, and accepts an injected helper so a test can stub both sides identically.

- [ ] **Step 1: Write the failing test**

```groovy
// modules/nextflow/src/test/groovy/nextflow/processor/hash/HashContextTest.groovy
package nextflow.processor.hash

import java.nio.file.Paths

import nextflow.Session
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

class HashContextTest extends Specification {

    def 'exposes task, processor and session, and delegates the helpers'() {
        given:
        def session = Mock(Session)
        def processor = Mock(TaskProcessor) {
            getSession() >> session
        }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
            getSource() >> 'echo hello'
        }
        and:
        def helper = Mock(TaskHasher) {
            getTaskGlobalVars() >> [foo: 'a']
            getTaskBinEntries('echo hello') >> [Paths.get('/bin/x.sh')]
        }
        and:
        def ctx = new HashContext(task, helper)

        expect:
        ctx.task.is(task)
        ctx.processor.is(processor)
        ctx.session.is(session)
        ctx.globalVars() == [foo: 'a']
        ctx.binEntries() == [Paths.get('/bin/x.sh')]
    }

    def 'the vocabulary has exactly the keys master hashes'() {
        expect:
        HashKey.values()*.name() as Set == [
            'SESSION_ID', 'PROCESS_NAME', 'TASK_SOURCE', 'CONTAINER', 'INPUTS',
            'EVAL_OUTPUTS', 'SCRIPT_VARS', 'BIN_ENTRIES', 'MODULE_BUNDLE',
            'ENV_MODULES', 'CONDA', 'SPACK', 'STUB_MARKER'
        ] as Set
    }
}
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.HashContextTest"`
Expected: FAIL — `unable to resolve class HashContext`.

- [ ] **Step 3: Create the four types**

```groovy
// modules/nextflow/src/main/groovy/nextflow/processor/hash/HashKey.groovy
package nextflow.processor.hash

import groovy.transform.CompileStatic

/**
 * The canonical vocabulary of task hash keys.
 *
 * Names are load-bearing data: a consumer comparing two runs aligns them by key
 * name, so a rename is a semantic change and must move the spec fingerprint.
 * The set is closed — a new key requires a constant here, which is what makes
 * the exhaustiveness check in TaskHashSpecTest possible.
 */
@CompileStatic
enum HashKey {
    SESSION_ID,
    PROCESS_NAME,
    TASK_SOURCE,
    CONTAINER,
    INPUTS,
    EVAL_OUTPUTS,
    SCRIPT_VARS,
    BIN_ENTRIES,
    MODULE_BUNDLE,
    ENV_MODULES,
    CONDA,
    SPACK,
    STUB_MARKER
}
```

```groovy
// modules/nextflow/src/main/groovy/nextflow/processor/hash/Contributor.groovy
package nextflow.processor.hash

/**
 * Extracts the value(s) one hash key contributes for a task.
 *
 * A contributor chooses WHICH objects enter the hash; EncodingRules decide HOW any
 * object becomes bytes. Emitting the wrong NUMBER of values breaks byte-exactness
 * just as surely as emitting the wrong value, so an absent key must emit an empty
 * list — never a list containing null.
 */
interface Contributor {

    /**
     * Stable identity of this extraction, contributing to the spec fingerprint.
     * Two specs that extract differently for the same key must differ here.
     */
    String canonicalName()

    List<Object> emit(HashContext ctx)
}
```

```groovy
// modules/nextflow/src/main/groovy/nextflow/processor/hash/KeyBinding.groovy
package nextflow.processor.hash

import groovy.transform.Canonical
import groovy.transform.CompileStatic

@Canonical
@CompileStatic
class KeyBinding {
    HashKey key
    Contributor contributor
}
```

```groovy
// modules/nextflow/src/main/groovy/nextflow/processor/hash/HashContext.groovy
package nextflow.processor.hash

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun

/**
 * Everything a Contributor may read while extracting values for a task.
 *
 * The two derived helpers are delegated to the existing TaskHasher rather than
 * duplicated, so the interpreter and the byte-exactness oracle cannot diverge on
 * them. The helper is injectable so a test can stub both sides identically.
 */
@CompileStatic
class HashContext {

    final TaskRun task

    final TaskProcessor processor

    final Session session

    /** Services injected by a plugin-supplied spec, keyed by name. */
    final Map<String,Object> services

    private final TaskHasher helper

    HashContext(TaskRun task) {
        this(task, new TaskHasher(task), [:])
    }

    HashContext(TaskRun task, TaskHasher helper) {
        this(task, helper, [:])
    }

    HashContext(TaskRun task, TaskHasher helper, Map<String,Object> services) {
        this.task = task
        this.processor = task.processor
        this.session = task.processor.session
        this.helper = helper
        this.services = services
    }

    Map<String,Object> globalVars() {
        return helper.getTaskGlobalVars()
    }

    List<Path> binEntries() {
        return helper.getTaskBinEntries(task.source)
    }
}
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.HashContextTest"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/hash/ \
        modules/nextflow/src/test/groovy/nextflow/processor/hash/
git commit -s -m "Add task hash key vocabulary and contributor types"
```

---

### Task 3: TaskHashSpec and its fingerprint

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/hash/TaskHashSpec.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/processor/hash/TaskHashSpecTest.groovy`

**Interfaces:**
- Consumes: `HashKey`, `KeyBinding`, `Contributor` (Task 2); `EncodingRules` (Task 1).
- Produces: `TaskHashSpec(String id, List<KeyBinding> bindings, EncodingRules encoding)` with `getId()`, `getBindings()`, `getEncoding()`, `canonicalForm()` returning `String`, `fingerprint()` returning `String`, and `keys()` returning `List<HashKey>`.

- [ ] **Step 1: Write the failing test**

```groovy
// modules/nextflow/src/test/groovy/nextflow/processor/hash/TaskHashSpecTest.groovy
package nextflow.processor.hash

import nextflow.util.EncodingRules
import spock.lang.Specification

class TaskHashSpecTest extends Specification {

    static Contributor contrib(String name) {
        return new Contributor() {
            @Override
            String canonicalName() { return name }
            @Override
            List<Object> emit(HashContext ctx) { return [] }
        }
    }

    def 'canonical form names the id, the ordered keys with their extraction, and the encoding'() {
        given:
        def spec = new TaskHashSpec('std/test', [
            new KeyBinding(HashKey.SESSION_ID, contrib('sessionId')),
            new KeyBinding(HashKey.TASK_SOURCE, contrib('taskSource'))
        ], EncodingRules.RECORD_TYPES)

        expect:
        spec.canonicalForm() == '''\
id=std/test
key=SESSION_ID:sessionId
key=TASK_SOURCE:taskSource
encoding=orderIndependentMaps=true;cacheFunnelFirst=true
function=murmur3_128
'''
    }

    def 'fingerprint is stable, and moves when any part of the canonical form moves'() {
        given:
        def base = new TaskHashSpec('std/test',
            [new KeyBinding(HashKey.SESSION_ID, contrib('sessionId'))], EncodingRules.RECORD_TYPES)

        expect: 'stable across calls and across equal specs'
        base.fingerprint() == base.fingerprint()
        and:
        base.fingerprint() == new TaskHashSpec('std/test',
            [new KeyBinding(HashKey.SESSION_ID, contrib('sessionId'))], EncodingRules.RECORD_TYPES).fingerprint()

        and: 'a different extraction for the same key moves it'
        base.fingerprint() != new TaskHashSpec('std/test',
            [new KeyBinding(HashKey.SESSION_ID, contrib('sessionId.other'))], EncodingRules.RECORD_TYPES).fingerprint()

        and: 'different encoding moves it'
        base.fingerprint() != new TaskHashSpec('std/test',
            [new KeyBinding(HashKey.SESSION_ID, contrib('sessionId'))], EncodingRules.LEGACY).fingerprint()

        and: 'different key order moves it'
        new TaskHashSpec('std/test', [
            new KeyBinding(HashKey.SESSION_ID, contrib('a')),
            new KeyBinding(HashKey.TASK_SOURCE, contrib('b'))
        ], EncodingRules.RECORD_TYPES).fingerprint() != new TaskHashSpec('std/test', [
            new KeyBinding(HashKey.TASK_SOURCE, contrib('b')),
            new KeyBinding(HashKey.SESSION_ID, contrib('a'))
        ], EncodingRules.RECORD_TYPES).fingerprint()
    }

    def 'keys() reports the declared keys in order'() {
        given:
        def spec = new TaskHashSpec('std/test', [
            new KeyBinding(HashKey.TASK_SOURCE, contrib('a')),
            new KeyBinding(HashKey.CONDA, contrib('b'))
        ], EncodingRules.RECORD_TYPES)

        expect:
        spec.keys() == [HashKey.TASK_SOURCE, HashKey.CONDA]
    }
}
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.TaskHashSpecTest"`
Expected: FAIL — `unable to resolve class TaskHashSpec`.

- [ ] **Step 3: Implement TaskHashSpec**

```groovy
// modules/nextflow/src/main/groovy/nextflow/processor/hash/TaskHashSpec.groovy
package nextflow.processor.hash

import com.google.common.hash.Hashing
import groovy.transform.CompileStatic
import nextflow.util.EncodingRules

/**
 * A task hash version expressed as data: which keys are hashed, in what order, how
 * each is extracted, and how values are encoded.
 *
 * A spec is immutable and, once shipped, frozen: editing one changes the hashes of
 * every task already recorded under it. A new behaviour is a new spec.
 */
@CompileStatic
class TaskHashSpec {

    final String id

    final List<KeyBinding> bindings

    final EncodingRules encoding

    private volatile String fingerprintValue

    TaskHashSpec(String id, List<KeyBinding> bindings, EncodingRules encoding) {
        this.id = id
        this.bindings = Collections.unmodifiableList(new ArrayList<KeyBinding>(bindings))
        this.encoding = encoding
    }

    List<HashKey> keys() {
        return bindings.collect { KeyBinding it -> it.key }
    }

    /**
     * Stable text form from which the fingerprint is derived. It names only semantics
     * — id, ordered keys, extraction identity, encoding, hash function — so that
     * refactoring implementation detail cannot move a historical spec's identity.
     * Never reformat this method.
     */
    String canonicalForm() {
        final sb = new StringBuilder()
        sb.append('id=').append(id).append('\n')
        for( KeyBinding b : bindings ) {
            sb.append('key=').append(b.key.name()).append(':').append(b.contributor.canonicalName()).append('\n')
        }
        sb.append('encoding=').append(encoding.canonicalForm()).append('\n')
        sb.append('function=murmur3_128').append('\n')
        return sb.toString()
    }

    String fingerprint() {
        if( fingerprintValue == null ) {
            fingerprintValue = Hashing.murmur3_128().newHasher()
                .putUnencodedChars(canonicalForm())
                .hash()
                .toString()
        }
        return fingerprintValue
    }

    @Override
    String toString() {
        return "TaskHashSpec[${id}@${fingerprint()}]"
    }
}
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.TaskHashSpecTest"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/hash/TaskHashSpec.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/hash/TaskHashSpecTest.groovy
git commit -s -m "Add TaskHashSpec with fingerprint derived from its canonical form"
```

---

### Task 4: The contributor library

Each contributor reproduces one key site of `TaskHasher.compute()` exactly, including its emission arity and its "emit nothing when absent" behaviour.

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/hash/Contributors.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/processor/hash/ContributorsTest.groovy`

**Interfaces:**
- Consumes: `HashContext`, `Contributor` (Task 2).
- Produces: static `Contributor` constants on `Contributors`: `SESSION_ID`, `PROCESS_NAME`, `TASK_SOURCE`, `CONTAINER`, `INPUTS_RAW`, `EVAL_OUTPUTS_RAW_MAP`, `EVAL_OUTPUTS_DERIVED_STRING`, `SCRIPT_VARS`, `BIN_ENTRIES`, `MODULE_BUNDLE`, `ENV_MODULES`, `CONDA`, `SPACK_AND_ARCH`, `STUB_MARKER`. Also `Contributors.of(String name, Closure<List<Object>> fn)` returning `Contributor`.

- [ ] **Step 1: Write the failing test**

```groovy
// modules/nextflow/src/test/groovy/nextflow/processor/hash/ContributorsTest.groovy
package nextflow.processor.hash

import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import spock.lang.Specification

class ContributorsTest extends Specification {

    private HashContext ctxFor(TaskRun task, TaskHasher helper) {
        return new HashContext(task, helper)
    }

    def 'absent optional keys emit nothing, not null'() {
        given:
        def session = Mock(Session)
        def processor = Mock(TaskProcessor) { getSession() >> session }
        def config = Mock(TaskConfig) {
            getModule() >> []
            getArchitecture() >> null
        }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
            isContainerEnabled() >> false
            getOutputEvals() >> [:]
            getCondaEnv() >> null
            getSpackEnv() >> null
            getConfig() >> config
        }
        def helper = Mock(TaskHasher)
        def ctx = ctxFor(task, helper)

        expect:
        Contributors.CONTAINER.emit(ctx) == []
        Contributors.EVAL_OUTPUTS_RAW_MAP.emit(ctx) == []
        Contributors.CONDA.emit(ctx) == []
        Contributors.SPACK_AND_ARCH.emit(ctx) == []
        Contributors.ENV_MODULES.emit(ctx) == []
    }

    def 'spack emits arch only when spack itself is set'() {
        given:
        def session = Mock(Session)
        def processor = Mock(TaskProcessor) { getSession() >> session }
        def config = Mock(TaskConfig) { getArchitecture() >> 'linux/amd64' }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
            getSpackEnv() >> spack
            getConfig() >> config
        }
        def ctx = ctxFor(task, Mock(TaskHasher))

        expect:
        Contributors.SPACK_AND_ARCH.emit(ctx) == expected

        where:
        spack       | expected
        null        | []
        'env.yaml'  | ['env.yaml', 'linux/amd64']
    }

    def 'inputs emit a name and a value per input, in order'() {
        given:
        def session = Mock(Session)
        def processor = Mock(TaskProcessor) { getSession() >> session }
        def task = Mock(TaskRun) {
            getProcessor() >> processor
            getInputs() >> [
                (Mock(nextflow.script.params.InParam) { getName() >> 'x' }): 1,
                (Mock(nextflow.script.params.InParam) { getName() >> 'y' }): 2
            ]
        }
        def ctx = ctxFor(task, Mock(TaskHasher))

        expect:
        Contributors.INPUTS_RAW.emit(ctx) == ['x', 1, 'y', 2]
    }

    def 'the two eval forms differ, and the derived form is the pre-7575 string'() {
        given:
        def session = Mock(Session)
        def processor = Mock(TaskProcessor) { getSession() >> session }
        def evals = [beta: 'echo b', alpha: 'echo a']
        def task = Mock(TaskRun) {
            getProcessor() >> processor
            getOutputEvals() >> evals
        }
        def ctx = ctxFor(task, Mock(TaskHasher))

        expect:
        Contributors.EVAL_OUTPUTS_RAW_MAP.emit(ctx) == ['eval_outputs', evals]
        and:
        Contributors.EVAL_OUTPUTS_DERIVED_STRING.emit(ctx) == ['eval_outputs', 'alpha=echo a\nbeta=echo b']
    }

    def 'canonical names are distinct across all contributors'() {
        given:
        def all = [
            Contributors.SESSION_ID, Contributors.PROCESS_NAME, Contributors.TASK_SOURCE,
            Contributors.CONTAINER, Contributors.INPUTS_RAW, Contributors.EVAL_OUTPUTS_RAW_MAP,
            Contributors.EVAL_OUTPUTS_DERIVED_STRING, Contributors.SCRIPT_VARS,
            Contributors.BIN_ENTRIES, Contributors.MODULE_BUNDLE, Contributors.ENV_MODULES,
            Contributors.CONDA, Contributors.SPACK_AND_ARCH, Contributors.STUB_MARKER
        ]

        expect:
        all*.canonicalName().toSet().size() == all.size()
    }
}
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.ContributorsTest"`
Expected: FAIL — `unable to resolve class Contributors`.

- [ ] **Step 3: Implement Contributors**

Each body is copied from the matching site in `TaskHasher.compute()`; do not "improve" any of them.

```groovy
// modules/nextflow/src/main/groovy/nextflow/processor/hash/Contributors.groovy
package nextflow.processor.hash

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.script.bundle.ResourcesBundle

/**
 * The contributor implementations, one per key site of TaskHasher.compute().
 *
 * Every body here mirrors the corresponding lines of that method exactly, including
 * how many values it appends and the fact that it appends none when the key is
 * absent. Byte-exactness lives in this file.
 */
@CompileStatic
class Contributors {

    static Contributor of(String name, Closure<List<Object>> fn) {
        return new Contributor() {
            @Override
            String canonicalName() {
                return name
            }

            @Override
            List<Object> emit(HashContext ctx) {
                return (List<Object>) fn.call(ctx)
            }
        }
    }

    static final Contributor SESSION_ID = of('sessionId') { HashContext ctx ->
        return [ctx.session.uniqueId] as List<Object>
    }

    static final Contributor PROCESS_NAME = of('processName') { HashContext ctx ->
        return [ctx.processor.name] as List<Object>
    }

    static final Contributor TASK_SOURCE = of('taskSource') { HashContext ctx ->
        return [ctx.task.source] as List<Object>
    }

    static final Contributor CONTAINER = of('containerFingerprint') { HashContext ctx ->
        if( !ctx.task.isContainerEnabled() ) {
            return [] as List<Object>
        }
        return [ctx.task.getContainerFingerprint()] as List<Object>
    }

    static final Contributor INPUTS_RAW = of('inputs.raw') { HashContext ctx ->
        final out = new ArrayList<Object>()
        for( final entry : ctx.task.inputs ) {
            out.add(entry.key.name)
            out.add(entry.value)
        }
        return out
    }

    /** Post-#7575: the eval map is hashed directly. */
    static final Contributor EVAL_OUTPUTS_RAW_MAP = of('evalOutputs.rawMap') { HashContext ctx ->
        final outEvals = ctx.task.getOutputEvals()
        if( !outEvals ) {
            return [] as List<Object>
        }
        return ['eval_outputs', outEvals] as List<Object>
    }

    /** Pre-#7575: a sorted "name=command" string, one entry per line. */
    static final Contributor EVAL_OUTPUTS_DERIVED_STRING = of('evalOutputs.derivedString') { HashContext ctx ->
        final outEvals = ctx.task.getOutputEvals()
        if( !outEvals ) {
            return [] as List<Object>
        }
        return ['eval_outputs', derivedEvalCommands(outEvals)] as List<Object>
    }

    static final Contributor SCRIPT_VARS = of('scriptVars') { HashContext ctx ->
        final vars = ctx.globalVars()
        if( !vars ) {
            return [] as List<Object>
        }
        return [vars.entrySet()] as List<Object>
    }

    static final Contributor BIN_ENTRIES = of('binEntries') { HashContext ctx ->
        final entries = ctx.binEntries()
        if( !entries ) {
            return [] as List<Object>
        }
        return new ArrayList<Object>(entries)
    }

    static final Contributor MODULE_BUNDLE = of('moduleBundleFingerprint') { HashContext ctx ->
        final ResourcesBundle bundle = ctx.session.enableModuleBinaries()
            ? ctx.processor.getModuleBundle()
            : null
        if( !bundle || !bundle.hasEntries() ) {
            return [] as List<Object>
        }
        return [bundle.fingerprint()] as List<Object>
    }

    static final Contributor ENV_MODULES = of('envModules') { HashContext ctx ->
        final modules = ctx.task.getConfig().getModule()
        if( !modules ) {
            return [] as List<Object>
        }
        return new ArrayList<Object>(modules)
    }

    static final Contributor CONDA = of('condaEnv') { HashContext ctx ->
        final conda = ctx.task.getCondaEnv()
        if( !conda ) {
            return [] as List<Object>
        }
        return [conda] as List<Object>
    }

    /** arch contributes only when spack is set — it is nested inside that branch today. */
    static final Contributor SPACK_AND_ARCH = of('spackEnvAndArch') { HashContext ctx ->
        final spack = ctx.task.getSpackEnv()
        if( !spack ) {
            return [] as List<Object>
        }
        final arch = ctx.task.getConfig().getArchitecture()
        if( !arch ) {
            return [spack] as List<Object>
        }
        return [spack, arch] as List<Object>
    }

    static final Contributor STUB_MARKER = of('stubMarker') { HashContext ctx ->
        if( ctx.session.stubRun && ctx.task.config.getStubBlock() ) {
            return ['stub-run'] as List<Object>
        }
        return [] as List<Object>
    }

    /**
     * Reproduces TaskHasher.computeEvalOutputCommands(), removed by #7575.
     * Retained verbatim because specs std/v1..v3 depend on its exact output.
     */
    protected static String derivedEvalCommands(Map<String,String> outEvals) {
        final result = new StringBuilder()
        final sortedEntries = outEvals.entrySet().sort { a, b -> a.key.compareTo(b.key) }
        for( final entry : sortedEntries ) {
            if( result.length() > 0 ) {
                result.append('\n')
            }
            result.append(entry.key).append('=').append(entry.value)
        }
        return result.toString()
    }
}
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.ContributorsTest"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/hash/Contributors.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/hash/ContributorsTest.groovy
git commit -s -m "Add task hash contributors mirroring TaskHasher key sites"
```

---

### Task 5: BaseTaskHasher, std/v4, and the byte-exactness oracle

This is the linchpin: `std/v4` must reproduce master's `TaskHasher` exactly.

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/hash/BaseTaskHasher.groovy`
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/hash/StdSpecs.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/processor/hash/BaseTaskHasherTest.groovy`

**Interfaces:**
- Consumes: `TaskHashSpec` (Task 3), `Contributors` (Task 4), `EncodingRules` (Task 1).
- Produces: `BaseTaskHasher(HashContext ctx, TaskHashSpec spec)` with `compute()` returning `HashCode` and `collectKeys()` returning `List<Object>`. `StdSpecs.STD_V4`, `StdSpecs.byId(String)` returning `TaskHashSpec` or throwing `IllegalArgumentException`.

- [ ] **Step 1: Write the failing test**

```groovy
// modules/nextflow/src/test/groovy/nextflow/processor/hash/BaseTaskHasherTest.groovy
package nextflow.processor.hash

import java.nio.file.Paths

import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ProcessConfig
import spock.lang.Specification

class BaseTaskHasherTest extends Specification {

    /** A task exercising every key the unit level can reach. */
    private Map fixture() {
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> false
            getStubRun() >> false
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'PIPE:FOO'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig)
            getModuleBundle() >> null
        }
        def config = Mock(TaskConfig) {
            getModule() >> ['gcc/11']
            getArchitecture() >> 'linux/amd64'
            getStubBlock() >> null
            getHashMode() >> nextflow.util.CacheHelper.HashMode.STANDARD
        }
        def task = Mock(TaskRun) {
            getSource() >> 'echo hello'
            getProcessor() >> processor
            getConfig() >> config
            isContainerEnabled() >> false
            getInputs() >> [:]
            getOutputEvals() >> [alpha: 'echo a']
            getCondaEnv() >> 'env.yml'
            getSpackEnv() >> 'spack.yaml'
        }
        def helper = Spy(new TaskHasher(task))
        helper.getTaskGlobalVars() >> [foo: 'a', bar: 'b']
        helper.getTaskBinEntries(_) >> []
        return [task: task, helper: helper]
    }

    def 'std/v4 reproduces the legacy TaskHasher byte for byte'() {
        given:
        def f = fixture()
        def legacy = f.helper as TaskHasher
        def ctx = new HashContext(f.task as TaskRun, legacy)

        when:
        def expected = legacy.compute()
        def actual = new BaseTaskHasher(ctx, StdSpecs.STD_V4).compute()

        then:
        actual == expected
    }

    def 'collectKeys reproduces the legacy key list, element for element'() {
        given:
        def f = fixture()
        def ctx = new HashContext(f.task as TaskRun, f.helper as TaskHasher)

        expect:
        new BaseTaskHasher(ctx, StdSpecs.STD_V4).collectKeys() == [
            UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c'),
            'PIPE:FOO',
            'echo hello',
            'eval_outputs', [alpha: 'echo a'],
            [foo: 'a', bar: 'b'].entrySet(),
            'gcc/11',
            'env.yml',
            'spack.yaml', 'linux/amd64'
        ]
    }

    def 'byId resolves known specs and rejects unknown ones'() {
        expect:
        StdSpecs.byId('std/v4').is(StdSpecs.STD_V4)

        when:
        StdSpecs.byId('std/nope')
        then:
        thrown(IllegalArgumentException)
    }
}
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.BaseTaskHasherTest"`
Expected: FAIL — `unable to resolve class BaseTaskHasher`.

- [ ] **Step 3: Implement BaseTaskHasher**

```groovy
// modules/nextflow/src/main/groovy/nextflow/processor/hash/BaseTaskHasher.groovy
package nextflow.processor.hash

import com.google.common.hash.HashCode
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.UnexpectedException
import nextflow.util.CacheHelper
import nextflow.util.HashBuilder

/**
 * Executes any TaskHashSpec. There are no subclasses: a different hash version, or
 * a plugin's own hasher, is a different spec value.
 */
@Slf4j
@CompileStatic
class BaseTaskHasher {

    private final HashContext ctx

    private final TaskHashSpec spec

    BaseTaskHasher(HashContext ctx, TaskHashSpec spec) {
        this.ctx = ctx
        this.spec = spec
    }

    TaskHashSpec getSpec() {
        return spec
    }

    /** The flat, ordered value list the spec produces for this task. */
    List<Object> collectKeys() {
        final keys = new ArrayList<Object>()
        for( KeyBinding binding : spec.bindings ) {
            keys.addAll(binding.contributor.emit(ctx))
        }
        return keys
    }

    HashCode compute() {
        final keys = collectKeys()
        final mode = ctx.task.processor.getConfig().getHashMode()
        try {
            return spec.encoding
                .apply(new HashBuilder().withHasher(HashBuilder.defaultHasher()).withMode(mode))
                .with(keys)
                .build()
        }
        catch( Throwable e ) {
            final msg = "Something went wrong while creating task hash for process '${ctx.processor.name}' under spec '${spec.id}' -- Offending keys: ${ keys.collect { k -> "\n - type=${k?.getClass()?.getName()} value=$k" } }"
            throw new UnexpectedException(msg, e)
        }
    }
}
```

Note: `CacheHelper.hasher(keys, mode)` is *not* used, because it offers no way to apply `EncodingRules`. `HashBuilder.with(List)` walks a `Collection` element by element, which is what `CacheHelper.hasher` does for the key list today — the `collectKeys` test in Step 1 plus the `std/v4` oracle test together confirm the two paths agree.

- [ ] **Step 4: Implement StdSpecs with std/v4 only**

```groovy
// modules/nextflow/src/main/groovy/nextflow/processor/hash/StdSpecs.groovy
package nextflow.processor.hash

import groovy.transform.CompileStatic
import nextflow.util.EncodingRules

/**
 * The standard task hash specs, one per historical behaviour of master's hash.
 *
 * These are frozen. Changing a published spec changes the hash of every task already
 * recorded under it; a new behaviour gets a new id.
 */
@CompileStatic
class StdSpecs {

    /** Current master: post-#7575 (2026-09-03), eval hashed as a raw map. */
    static final TaskHashSpec STD_V4 = new TaskHashSpec('std/v4', [
        new KeyBinding(HashKey.SESSION_ID, Contributors.SESSION_ID),
        new KeyBinding(HashKey.PROCESS_NAME, Contributors.PROCESS_NAME),
        new KeyBinding(HashKey.TASK_SOURCE, Contributors.TASK_SOURCE),
        new KeyBinding(HashKey.CONTAINER, Contributors.CONTAINER),
        new KeyBinding(HashKey.INPUTS, Contributors.INPUTS_RAW),
        new KeyBinding(HashKey.EVAL_OUTPUTS, Contributors.EVAL_OUTPUTS_RAW_MAP),
        new KeyBinding(HashKey.SCRIPT_VARS, Contributors.SCRIPT_VARS),
        new KeyBinding(HashKey.BIN_ENTRIES, Contributors.BIN_ENTRIES),
        new KeyBinding(HashKey.MODULE_BUNDLE, Contributors.MODULE_BUNDLE),
        new KeyBinding(HashKey.ENV_MODULES, Contributors.ENV_MODULES),
        new KeyBinding(HashKey.CONDA, Contributors.CONDA),
        new KeyBinding(HashKey.SPACK, Contributors.SPACK_AND_ARCH),
        new KeyBinding(HashKey.STUB_MARKER, Contributors.STUB_MARKER)
    ], EncodingRules.RECORD_TYPES)

    static final TaskHashSpec DEFAULT = STD_V4

    static List<TaskHashSpec> all() {
        return [STD_V4]
    }

    static TaskHashSpec byId(String id) {
        final found = all().find { TaskHashSpec it -> it.id == id }
        if( !found ) {
            throw new IllegalArgumentException("Unknown task hash spec: ${id} -- available: ${all()*.id.join(', ')}")
        }
        return found
    }
}
```

- [ ] **Step 5: Run the test and confirm it passes**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.BaseTaskHasherTest"`
Expected: PASS — in particular `std/v4 reproduces the legacy TaskHasher byte for byte`.

If the oracle test fails, do not adjust `StdSpecs`. Print both key lists (`collectKeys()` against a debug dump of `TaskHasher.compute()`'s local `keys`) and fix the contributor whose emission differs — the fault is always arity or a conditional, never the fold.

- [ ] **Step 6: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/hash/BaseTaskHasher.groovy \
        modules/nextflow/src/main/groovy/nextflow/processor/hash/StdSpecs.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/hash/BaseTaskHasherTest.groovy
git commit -s -m "Add BaseTaskHasher interpreter and the std/v4 spec"
```

---

### Task 6: The historical specs std/v1, std/v2, std/v3

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/processor/hash/StdSpecs.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/processor/hash/StdSpecsTest.groovy`

**Interfaces:**
- Consumes: everything from Task 5.
- Produces: `StdSpecs.STD_V1`, `StdSpecs.STD_V2`, `StdSpecs.STD_V3`; `StdSpecs.all()` returns all four in order.

Composition, from the verified history in the spec:

| | `MODULE_BUNDLE` | `EVAL_OUTPUTS` | encoding |
| --- | --- | --- | --- |
| `std/v1` | absent | derived string | `LEGACY` |
| `std/v2` | absent | derived string | `RECORD_TYPES` |
| `std/v3` | present | derived string | `RECORD_TYPES` |
| `std/v4` | present | raw map | `RECORD_TYPES` |

- [ ] **Step 1: Write the failing test**

```groovy
// modules/nextflow/src/test/groovy/nextflow/processor/hash/StdSpecsTest.groovy
package nextflow.processor.hash

import nextflow.util.EncodingRules
import spock.lang.Specification

class StdSpecsTest extends Specification {

    def 'adjacent specs differ in exactly the documented place'() {
        expect: 'v1 to v2 differs only in encoding'
        StdSpecs.STD_V1.keys() == StdSpecs.STD_V2.keys()
        StdSpecs.STD_V1.encoding.is(EncodingRules.LEGACY)
        StdSpecs.STD_V2.encoding.is(EncodingRules.RECORD_TYPES)

        and: 'v2 to v3 differs only by the module bundle key'
        StdSpecs.STD_V3.keys() - StdSpecs.STD_V2.keys() == [HashKey.MODULE_BUNDLE]
        StdSpecs.STD_V2.keys() - StdSpecs.STD_V3.keys() == []
        StdSpecs.STD_V3.encoding.is(StdSpecs.STD_V2.encoding)

        and: 'v3 to v4 differs only in the eval extraction'
        StdSpecs.STD_V3.keys() == StdSpecs.STD_V4.keys()
        StdSpecs.STD_V3.bindings.find { it.key == HashKey.EVAL_OUTPUTS }.contributor.canonicalName() ==
            'evalOutputs.derivedString'
        StdSpecs.STD_V4.bindings.find { it.key == HashKey.EVAL_OUTPUTS }.contributor.canonicalName() ==
            'evalOutputs.rawMap'
    }

    def 'bin entries are present in every spec, as they have been since 2018'() {
        expect:
        StdSpecs.all().every { HashKey.BIN_ENTRIES in it.keys() }
    }

    def 'every spec accounts for every key in the vocabulary'() {
        given: 'keys a spec may legitimately omit, with the reason'
        def permittedOmissions = [
            'std/v1': [HashKey.MODULE_BUNDLE] as Set,   // added by #6914, 2026-07-17
            'std/v2': [HashKey.MODULE_BUNDLE] as Set,
            'std/v3': [] as Set,
            'std/v4': [] as Set
        ]

        expect:
        StdSpecs.all().every { spec ->
            def missing = (HashKey.values() as Set) - (spec.keys() as Set)
            missing == permittedOmissions[spec.id]
        }
    }

    def 'fingerprints are distinct across all specs'() {
        expect:
        StdSpecs.all()*.fingerprint().toSet().size() == StdSpecs.all().size()
    }
}
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.StdSpecsTest"`
Expected: FAIL — `No such property: STD_V1`.

- [ ] **Step 3: Add the three historical specs**

Insert before `STD_V4` in `StdSpecs.groovy`:

```groovy
    /** Before #6679 "Record types" (2026-03-09): legacy encoding, no module bundle key. */
    static final TaskHashSpec STD_V1 = new TaskHashSpec('std/v1', [
        new KeyBinding(HashKey.SESSION_ID, Contributors.SESSION_ID),
        new KeyBinding(HashKey.PROCESS_NAME, Contributors.PROCESS_NAME),
        new KeyBinding(HashKey.TASK_SOURCE, Contributors.TASK_SOURCE),
        new KeyBinding(HashKey.CONTAINER, Contributors.CONTAINER),
        new KeyBinding(HashKey.INPUTS, Contributors.INPUTS_RAW),
        new KeyBinding(HashKey.EVAL_OUTPUTS, Contributors.EVAL_OUTPUTS_DERIVED_STRING),
        new KeyBinding(HashKey.SCRIPT_VARS, Contributors.SCRIPT_VARS),
        new KeyBinding(HashKey.BIN_ENTRIES, Contributors.BIN_ENTRIES),
        new KeyBinding(HashKey.ENV_MODULES, Contributors.ENV_MODULES),
        new KeyBinding(HashKey.CONDA, Contributors.CONDA),
        new KeyBinding(HashKey.SPACK, Contributors.SPACK_AND_ARCH),
        new KeyBinding(HashKey.STUB_MARKER, Contributors.STUB_MARKER)
    ], EncodingRules.LEGACY)

    /** #6679 (2026-03-09) to #6914 (2026-07-17): record-types encoding, still no module bundle. */
    static final TaskHashSpec STD_V2 = new TaskHashSpec('std/v2', STD_V1.bindings, EncodingRules.RECORD_TYPES)

    /** #6914 (2026-07-17) to #7575 (2026-09-03): module bundle key added. */
    static final TaskHashSpec STD_V3 = new TaskHashSpec('std/v3', [
        new KeyBinding(HashKey.SESSION_ID, Contributors.SESSION_ID),
        new KeyBinding(HashKey.PROCESS_NAME, Contributors.PROCESS_NAME),
        new KeyBinding(HashKey.TASK_SOURCE, Contributors.TASK_SOURCE),
        new KeyBinding(HashKey.CONTAINER, Contributors.CONTAINER),
        new KeyBinding(HashKey.INPUTS, Contributors.INPUTS_RAW),
        new KeyBinding(HashKey.EVAL_OUTPUTS, Contributors.EVAL_OUTPUTS_DERIVED_STRING),
        new KeyBinding(HashKey.SCRIPT_VARS, Contributors.SCRIPT_VARS),
        new KeyBinding(HashKey.BIN_ENTRIES, Contributors.BIN_ENTRIES),
        new KeyBinding(HashKey.MODULE_BUNDLE, Contributors.MODULE_BUNDLE),
        new KeyBinding(HashKey.ENV_MODULES, Contributors.ENV_MODULES),
        new KeyBinding(HashKey.CONDA, Contributors.CONDA),
        new KeyBinding(HashKey.SPACK, Contributors.SPACK_AND_ARCH),
        new KeyBinding(HashKey.STUB_MARKER, Contributors.STUB_MARKER)
    ], EncodingRules.RECORD_TYPES)
```

And replace `all()`:

```groovy
    static List<TaskHashSpec> all() {
        return [STD_V1, STD_V2, STD_V3, STD_V4]
    }
```

- [ ] **Step 4: Run the tests and confirm they pass**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.*"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/hash/StdSpecs.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/hash/StdSpecsTest.groovy
git commit -s -m "Add historical task hash specs std/v1 to std/v3"
```

---

### Task 7: Per-key digests

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/processor/hash/BaseTaskHasher.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/processor/hash/BaseTaskHasherExplainTest.groovy`

**Interfaces:**
- Consumes: `BaseTaskHasher` (Task 5).
- Produces: `BaseTaskHasher.explain()` returning `Map<HashKey,HashCode>` — insertion-ordered, one entry per binding that emitted at least one value.

- [ ] **Step 1: Write the failing test**

```groovy
// modules/nextflow/src/test/groovy/nextflow/processor/hash/BaseTaskHasherExplainTest.groovy
package nextflow.processor.hash

import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ProcessConfig
import spock.lang.Specification

class BaseTaskHasherExplainTest extends Specification {

    private HashContext ctx(String source, String conda) {
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> false
            getStubRun() >> false
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'PIPE:FOO'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig)
            getModuleBundle() >> null
        }
        def config = Mock(TaskConfig) {
            getModule() >> []
            getArchitecture() >> null
            getStubBlock() >> null
            getHashMode() >> nextflow.util.CacheHelper.HashMode.STANDARD
        }
        def task = Mock(TaskRun) {
            getSource() >> source
            getProcessor() >> processor
            getConfig() >> config
            isContainerEnabled() >> false
            getInputs() >> [:]
            getOutputEvals() >> [:]
            getCondaEnv() >> conda
            getSpackEnv() >> null
        }
        def helper = Spy(new TaskHasher(task))
        helper.getTaskGlobalVars() >> [:]
        helper.getTaskBinEntries(_) >> []
        return new HashContext(task, helper)
    }

    def 'explain reports one digest per emitting key, in spec order'() {
        when:
        def explained = new BaseTaskHasher(ctx('echo a', 'env.yml'), StdSpecs.STD_V4).explain()

        then: 'keys that emitted nothing are absent'
        explained.keySet() as List == [
            HashKey.SESSION_ID, HashKey.PROCESS_NAME, HashKey.TASK_SOURCE, HashKey.CONDA
        ]
    }

    def 'only the differing key changes digest between two tasks'() {
        given:
        def a = new BaseTaskHasher(ctx('echo a', 'env.yml'), StdSpecs.STD_V4).explain()
        def b = new BaseTaskHasher(ctx('echo CHANGED', 'env.yml'), StdSpecs.STD_V4).explain()

        expect:
        a[HashKey.TASK_SOURCE] != b[HashKey.TASK_SOURCE]
        and:
        (a.keySet() - HashKey.TASK_SOURCE).every { a[it] == b[it] }
    }

    def 'explain does not change the computed hash'() {
        given:
        def hasher = new BaseTaskHasher(ctx('echo a', 'env.yml'), StdSpecs.STD_V4)
        def before = hasher.compute()

        when:
        hasher.explain()

        then:
        hasher.compute() == before
    }
}
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.BaseTaskHasherExplainTest"`
Expected: FAIL — `No signature of method: explain()`.

- [ ] **Step 3: Add explain()**

Append to `BaseTaskHasher`:

```groovy
    /**
     * Per-key digest of this task under this spec, for explaining a cache miss.
     *
     * Strictly additive: it recomputes nothing that feeds compute(), and a key that
     * emits no value is omitted rather than digested as empty.
     */
    Map<HashKey,HashCode> explain() {
        final mode = ctx.task.processor.getConfig().getHashMode()
        final result = new LinkedHashMap<HashKey,HashCode>()
        for( KeyBinding binding : spec.bindings ) {
            final values = binding.contributor.emit(ctx)
            if( !values ) {
                continue
            }
            result.put(binding.key, spec.encoding
                .apply(new HashBuilder().withHasher(HashBuilder.defaultHasher()).withMode(mode))
                .with(values)
                .build())
        }
        return result
    }
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.BaseTaskHasherExplainTest"`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/hash/BaseTaskHasher.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/hash/BaseTaskHasherExplainTest.groovy
git commit -s -m "Add per-key digests to BaseTaskHasher"
```

---

### Task 8: Spec selection and TaskProcessor wiring

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/hash/TaskHashSpecFactory.groovy`
- Create: `modules/nextflow/src/main/groovy/nextflow/processor/hash/TaskHashSpecResolver.groovy`
- Modify: `modules/nextflow/src/main/groovy/nextflow/processor/TaskProcessor.groovy:684`
- Test: `modules/nextflow/src/test/groovy/nextflow/processor/hash/TaskHashSpecResolverTest.groovy`

**Interfaces:**
- Consumes: `StdSpecs` (Tasks 5–6), `TaskHashSpec` (Task 3).
- Produces: `interface TaskHashSpecFactory extends ExtensionPoint { TaskHashSpec create(TaskRun task) }`; `TaskHashSpecResolver.resolve(TaskRun task)` returning `TaskHashSpec`; `TaskHashSpecResolver.defaultSpec()` returning `TaskHashSpec`.

Resolution order: plugin factories in priority order (first non-null wins), then `NXF_TASK_HASH_VER` if set, then `StdSpecs.DEFAULT`.

- [ ] **Step 1: Write the failing test**

```groovy
// modules/nextflow/src/test/groovy/nextflow/processor/hash/TaskHashSpecResolverTest.groovy
package nextflow.processor.hash

import nextflow.SysEnv
import nextflow.processor.TaskRun
import spock.lang.Specification

class TaskHashSpecResolverTest extends Specification {

    def cleanup() {
        SysEnv.pop()
    }

    def 'defaults to std/v4 when nothing selects a spec'() {
        given:
        SysEnv.push([:])

        expect:
        TaskHashSpecResolver.defaultSpec().is(StdSpecs.STD_V4)
    }

    def 'NXF_TASK_HASH_VER selects a spec by id'() {
        given:
        SysEnv.push([NXF_TASK_HASH_VER: 'std/v2'])

        expect:
        TaskHashSpecResolver.defaultSpec().is(StdSpecs.STD_V2)
    }

    def 'an unknown id fails loudly rather than silently falling back'() {
        given:
        SysEnv.push([NXF_TASK_HASH_VER: 'std/nope'])

        when:
        TaskHashSpecResolver.defaultSpec()
        then:
        thrown(IllegalArgumentException)
    }

    def 'a plugin factory takes precedence over the default'() {
        given:
        SysEnv.push([:])
        def custom = new TaskHashSpec('plugin/v1', StdSpecs.STD_V4.bindings, StdSpecs.STD_V4.encoding)
        def factory = Mock(TaskHashSpecFactory)
        def task = Mock(TaskRun)

        when:
        def result = TaskHashSpecResolver.resolve(task, [factory])
        then:
        1 * factory.create(task) >> custom
        result.is(custom)
    }

    def 'a factory that abstains falls through to the next, then to the default'() {
        given:
        SysEnv.push([:])
        def abstaining = Mock(TaskHashSpecFactory)
        def task = Mock(TaskRun)

        when:
        def result = TaskHashSpecResolver.resolve(task, [abstaining])
        then:
        1 * abstaining.create(task) >> null
        result.is(StdSpecs.STD_V4)
    }
}
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.TaskHashSpecResolverTest"`
Expected: FAIL — `unable to resolve class TaskHashSpecFactory`.

- [ ] **Step 3: Create the extension point and resolver**

```groovy
// modules/nextflow/src/main/groovy/nextflow/processor/hash/TaskHashSpecFactory.groovy
package nextflow.processor.hash

import nextflow.processor.TaskRun
import org.pf4j.ExtensionPoint

/**
 * Extension point through which a plugin supplies the task hash spec for a task.
 *
 * Plugins contribute a spec rather than a hasher so that every hash goes through
 * BaseTaskHasher: per-key digests and spec fingerprints then hold for plugin-driven
 * runs too, and a plugin cannot duplicate the key list to vary two of its entries.
 */
interface TaskHashSpecFactory extends ExtensionPoint {

    /**
     * @return the spec for {@code task}, or {@code null} to abstain, in which case the
     *      next factory is asked and, failing all, the default spec is used.
     */
    TaskHashSpec create(TaskRun task)
}
```

```groovy
// modules/nextflow/src/main/groovy/nextflow/processor/hash/TaskHashSpecResolver.groovy
package nextflow.processor.hash

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import nextflow.SysEnv
import nextflow.plugin.Plugins
import nextflow.processor.TaskRun

/**
 * Chooses the spec for a task: plugin factories first, then NXF_TASK_HASH_VER, then
 * the current default.
 */
@CompileStatic
class TaskHashSpecResolver {

    static TaskHashSpec resolve(TaskRun task) {
        return resolve(task, Plugins.getPriorityExtensions(TaskHashSpecFactory) ?: Collections.<TaskHashSpecFactory>emptyList())
    }

    static TaskHashSpec resolve(TaskRun task, List<TaskHashSpecFactory> factories) {
        for( TaskHashSpecFactory factory : factories ) {
            final spec = factory.create(task)
            if( spec != null ) {
                return spec
            }
        }
        return defaultSpec()
    }

    static TaskHashSpec defaultSpec() {
        final id = SysEnv.get('NXF_TASK_HASH_VER')
        if( !id ) {
            return StdSpecs.DEFAULT
        }
        return StdSpecs.byId(id)
    }
}
```

- [ ] **Step 4: Run the test and confirm it passes**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.TaskHashSpecResolverTest"`
Expected: PASS.

- [ ] **Step 5: Wire it into TaskProcessor**

Replace line 684 of `modules/nextflow/src/main/groovy/nextflow/processor/TaskProcessor.groovy`:

```groovy
        final hash = new TaskHasher(task).compute()
```

with:

```groovy
        final spec = TaskHashSpecResolver.resolve(task)
        final hash = new BaseTaskHasher(new HashContext(task), spec).compute()
```

and add the imports beside the existing `nextflow.processor` imports:

```groovy
import nextflow.processor.hash.BaseTaskHasher
import nextflow.processor.hash.HashContext
import nextflow.processor.hash.TaskHashSpecResolver
```

- [ ] **Step 6: Run the full processor suite to confirm nothing moved**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.*"`
Expected: PASS. `TaskHasherTest` still passes because the class is untouched and remains the oracle.

- [ ] **Step 7: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/hash/TaskHashSpecFactory.groovy \
        modules/nextflow/src/main/groovy/nextflow/processor/hash/TaskHashSpecResolver.groovy \
        modules/nextflow/src/main/groovy/nextflow/processor/TaskProcessor.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/hash/TaskHashSpecResolverTest.groovy
git commit -s -m "Resolve the task hash spec per task and route hashing through it"
```

---

### Task 9: Named entries in `-dump-hashes json`

The immediate payoff: today's output has no key names, so the primary cache-miss debugging tool produces an unlabelled list.

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/processor/hash/BaseTaskHasher.groovy`
- Modify: `modules/nextflow/src/main/groovy/nextflow/processor/TaskProcessor.groovy` (the same site as Task 8)
- Test: `modules/nextflow/src/test/groovy/nextflow/processor/hash/BaseTaskHasherDumpTest.groovy`

**Interfaces:**
- Consumes: `BaseTaskHasher.explain()` (Task 7).
- Produces: `BaseTaskHasher.dumpJson()` returning `String` — a JSON array of `{key, hash}` objects plus a leading `{spec, fingerprint}` object.

- [ ] **Step 1: Write the failing test**

```groovy
// modules/nextflow/src/test/groovy/nextflow/processor/hash/BaseTaskHasherDumpTest.groovy
package nextflow.processor.hash

import groovy.json.JsonSlurper
import nextflow.Session
import nextflow.processor.TaskConfig
import nextflow.processor.TaskHasher
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.ProcessConfig
import spock.lang.Specification

class BaseTaskHasherDumpTest extends Specification {

    def 'dumpJson names every entry and identifies the spec'() {
        given:
        def session = Mock(Session) {
            getUniqueId() >> UUID.fromString('b69b6eeb-b332-4d2c-9957-c291b15f498c')
            enableModuleBinaries() >> false
            getStubRun() >> false
        }
        def processor = Mock(TaskProcessor) {
            getName() >> 'PIPE:FOO'
            getSession() >> session
            getConfig() >> Mock(ProcessConfig)
            getModuleBundle() >> null
        }
        def config = Mock(TaskConfig) {
            getModule() >> []
            getArchitecture() >> null
            getStubBlock() >> null
            getHashMode() >> nextflow.util.CacheHelper.HashMode.STANDARD
        }
        def task = Mock(TaskRun) {
            getSource() >> 'echo a'
            getProcessor() >> processor
            getConfig() >> config
            isContainerEnabled() >> false
            getInputs() >> [:]
            getOutputEvals() >> [:]
            getCondaEnv() >> null
            getSpackEnv() >> null
        }
        def helper = Spy(new TaskHasher(task))
        helper.getTaskGlobalVars() >> [:]
        helper.getTaskBinEntries(_) >> []

        when:
        def json = new JsonSlurper().parseText(
            new BaseTaskHasher(new HashContext(task, helper), StdSpecs.STD_V4).dumpJson()) as List

        then:
        json[0].spec == 'std/v4'
        json[0].fingerprint == StdSpecs.STD_V4.fingerprint()
        and:
        json[1..-1]*.key == ['SESSION_ID', 'PROCESS_NAME', 'TASK_SOURCE']
        and:
        json[1..-1].every { it.hash instanceof String && it.hash.length() > 0 }
    }
}
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.BaseTaskHasherDumpTest"`
Expected: FAIL — `No signature of method: dumpJson()`.

- [ ] **Step 3: Add dumpJson()**

Add the import `import groovy.json.JsonOutput` to `BaseTaskHasher`, then append:

```groovy
    /** Named per-key entries for `-dump-hashes json`, prefixed by the spec identity. */
    String dumpJson() {
        final entries = new ArrayList<Map<String,Object>>()
        entries.add([spec: spec.id, fingerprint: spec.fingerprint()] as Map<String,Object>)
        for( Map.Entry<HashKey,HashCode> e : explain().entrySet() ) {
            entries.add([key: e.key.name(), hash: e.value.toString()] as Map<String,Object>)
        }
        return JsonOutput.prettyPrint(JsonOutput.toJson(entries))
    }
```

- [ ] **Step 4: Emit it from the hashing site**

In `TaskProcessor.groovy`, extend the block added in Task 8:

```groovy
        final spec = TaskHashSpecResolver.resolve(task)
        final hasher = new BaseTaskHasher(new HashContext(task), spec)
        final hash = hasher.compute()
        if( session.dumpHashes == 'json' ) {
            log.info "[${task.name}] cache hash: ${hash}; spec: ${spec.id}; entries: ${hasher.dumpJson()}"
        }
```

- [ ] **Step 5: Run the tests and confirm they pass**

Run: `./gradlew :nextflow:test --tests "nextflow.processor.hash.*" --tests "nextflow.processor.TaskProcessorTest"`
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/processor/hash/BaseTaskHasher.groovy \
        modules/nextflow/src/main/groovy/nextflow/processor/TaskProcessor.groovy \
        modules/nextflow/src/test/groovy/nextflow/processor/hash/BaseTaskHasherDumpTest.groovy
git commit -s -m "Emit named per-key entries for -dump-hashes json"
```

---

### Task 10: Validation pipeline and runbook

The end-to-end proof: a cache hit on `-resume` against an older release *is* byte-exact hash equality.

**Files:**
- Create: `specs/260916-task-hash-spec/pipeline/main.nf`
- Create: `specs/260916-task-hash-spec/pipeline/nextflow.config`
- Create: `specs/260916-task-hash-spec/pipeline/bin/helper.sh`
- Create: `specs/260916-task-hash-spec/pipeline/modules/bundled/main.nf`
- Create: `specs/260916-task-hash-spec/pipeline/modules/bundled/resources/bin/tool.sh`
- Create: `specs/260916-task-hash-spec/pipeline/data/input.txt`
- Create: `specs/260916-task-hash-spec/runbook.md`

**Interfaces:**
- Consumes: `NXF_TASK_HASH_VER` (Task 8).
- Produces: nothing consumed by later tasks.

- [ ] **Step 1: Create the pipeline inputs**

```bash
mkdir -p specs/260916-task-hash-spec/pipeline/{bin,data,modules/bundled/resources/bin}
printf 'hello\n' > specs/260916-task-hash-spec/pipeline/data/input.txt
printf '#!/bin/bash\necho "helper says $1"\n' > specs/260916-task-hash-spec/pipeline/bin/helper.sh
printf '#!/bin/bash\necho "bundled tool"\n' > specs/260916-task-hash-spec/pipeline/modules/bundled/resources/bin/tool.sh
chmod +x specs/260916-task-hash-spec/pipeline/bin/helper.sh \
         specs/260916-task-hash-spec/pipeline/modules/bundled/resources/bin/tool.sh
```

- [ ] **Step 2: Write the pipeline**

```groovy
// specs/260916-task-hash-spec/pipeline/modules/bundled/main.nf
process P_MODULE_BUNDLE {
    output:
    stdout

    script:
    """
    tool.sh
    """
}
```

```groovy
// specs/260916-task-hash-spec/pipeline/main.nf
// One process per hash key, so a cache miss localises immediately.
// See ../runbook.md for how this is run.

include { P_MODULE_BUNDLE } from './modules/bundled/main.nf'

params.greeting = 'hello'

process P_BASIC {
    ext.flavour = 'vanilla'

    input:
    val word
    path sample

    output:
    stdout

    script:
    """
    echo "${word} ${params.greeting} ${task.ext.flavour} \$(cat ${sample})"
    """
}

// H1 vs H2 discriminator: a genuine Map reaches the hasher, so the encoding
// rules (values() vs entrySet()) actually bite.
process P_MAP_INPUT {
    input:
    val settings

    output:
    stdout

    script:
    """
    echo "${settings.alpha} ${settings.beta}"
    """
}

// H3 vs H4 discriminator: without an eval output the two specs are identical.
process P_EVAL {
    output:
    stdout
    eval 'echo evaluated', emit: probe

    script:
    """
    echo eval-process
    """
}

process P_BIN {
    output:
    stdout

    script:
    """
    helper.sh one
    """
}

process P_CONTAINER {
    container 'quay.io/nextflow/bash@sha256:ca9f6e1a3b1d2f3e4a5b6c7d8e9f0a1b2c3d4e5f60718293a4b5c6d7e8f90a1b'

    output:
    stdout

    script:
    """
    echo containerised
    """
}

process P_CONDA {
    conda 'bioconda::fastqc=0.12.1'

    output:
    stdout

    script:
    """
    echo conda-backed
    """
}

process P_STUB {
    output:
    stdout

    script:
    """
    echo real
    """

    stub:
    """
    echo stubbed
    """
}

workflow {
    P_BASIC(Channel.value('word'), file("${projectDir}/data/input.txt"))
    P_MAP_INPUT(Channel.value([alpha: 'one', beta: 'two']))
    P_EVAL()
    P_BIN()
    P_MODULE_BUNDLE()
    P_CONTAINER()
    P_CONDA()
    P_STUB()
}
```

```groovy
// specs/260916-task-hash-spec/pipeline/nextflow.config
// Module binaries must be on, or P_MODULE_BUNDLE never contributes the bundle
// fingerprint and the H2 vs H3 boundary goes untested.
nextflow.enable.moduleBinaries = true

profiles {
    withtooling {
        docker.enabled = true
        conda.enabled = true
    }
}
```

- [ ] **Step 3: Confirm the pipeline runs on master**

Run: `./launch.sh run specs/260916-task-hash-spec/pipeline/main.nf -profile withtooling`
Expected: all eight processes complete. If `P_CONTAINER` fails, replace the pinned digest with one that exists locally — but keep it pinned by digest, or the `CONTAINER` key varies for environmental reasons.

- [ ] **Step 4: Write the runbook**

````markdown
<!-- specs/260916-task-hash-spec/runbook.md -->
# Byte-exactness runbook

A cache hit on `-resume` is byte-exact hash equality, proven through `CacheDB` and
`TaskProcessor` rather than a unit assertion.

Run everything from `specs/260916-task-hash-spec/pipeline/`. `$PROTO` is the
prototype launcher (`../../../launch.sh`).

## 1. Positive: each spec reproduces its release

| Spec | Baseline |
| --- | --- |
| `std/v1` | `NXF_VER=26.02.0-edge` |
| `std/v2` | `NXF_VER=26.07.0-edge` |
| `std/v3` | `NXF_VER=26.08.0-edge` |
| `std/v4` | the prototype itself, no env var |

For each row:

```bash
rm -rf work .nextflow*
NXF_VER=<version> nextflow run main.nf -profile withtooling
NXF_TASK_HASH_VER=<spec> $PROTO run main.nf -profile withtooling -resume
```

**Pass:** every process reports `CACHED` in the second run.

## 2. Negative controls

Without these, a run where everything caches for unrelated reasons reads as a pass.
After a passing step 1, each of these must make exactly the named process re-execute:

```bash
# change a script body
sed -i 's/echo eval-process/echo eval-process-CHANGED/' main.nf   # P_EVAL re-runs
# change a value input
sed -i "s/Channel.value('word')/Channel.value('WORD')/" main.nf   # P_BASIC re-runs
# rename a process
sed -i 's/P_BIN/P_BIN_RENAMED/g' main.nf                          # P_BIN re-runs
```

Revert each change before the next.

## 3. Cross-spec discrimination

The strongest check: the specs must differ in exactly the predicted place.

```bash
rm -rf work .nextflow*
NXF_VER=26.08.0-edge nextflow run main.nf -profile withtooling      # H3
NXF_TASK_HASH_VER=std/v4 $PROTO run main.nf -profile withtooling -resume
```

**Pass:** everything `CACHED` **except `P_EVAL`**, which re-runs.

```bash
rm -rf work .nextflow*
NXF_VER=26.07.0-edge nextflow run main.nf -profile withtooling      # H2
NXF_TASK_HASH_VER=std/v3 $PROTO run main.nf -profile withtooling -resume
```

**Pass:** only `P_MODULE_BUNDLE` re-runs.

```bash
rm -rf work .nextflow*
NXF_VER=26.02.0-edge nextflow run main.nf -profile withtooling      # H1
NXF_TASK_HASH_VER=std/v2 $PROTO run main.nf -profile withtooling -resume
```

**Pass:** `P_MAP_INPUT` re-runs. If *nothing* re-runs, the Map never reached the
hasher and the encoding dimension is untested — fix the fixture before trusting
`std/v1`.

## 4. Stub marker

```bash
rm -rf work .nextflow*
NXF_VER=26.08.0-edge nextflow run main.nf -profile withtooling -stub-run
NXF_TASK_HASH_VER=std/v3 $PROTO run main.nf -profile withtooling -stub-run -resume
```

**Pass:** `P_STUB` reports `CACHED`.

## When resume fails

Resume can fail for reasons unrelated to hashing — a February build's LevelDB/Kryo
cache may not be readable by a master build. Use the independent second track to
tell the two apart:

```bash
NXF_VER=<version> nextflow run main.nf -profile withtooling -dump-hashes json 2>&1 | tee baseline.log
NXF_TASK_HASH_VER=<spec> $PROTO run main.nf -profile withtooling -dump-hashes json 2>&1 | tee proto.log
```

Compare the reported `cache hash:` per process. Equal hashes with a failing resume
means the cache plumbing, not the spec. Unequal hashes: the prototype's per-key
entries name the key that differs.

## Known limitation

The era map is derived from the history of `TaskHasher` and `HashBuilder` only. A
change in a `CacheFunnel` implementor (`FileHolder`, `ArrayBag`, …) would move hashes
without touching either file. If a positive run misses, either the spec is wrong or
there is a boundary we have not found — `-dump-hashes` says which key.
````

- [ ] **Step 5: Run step 1 of the runbook for `std/v4`**

Run: the `std/v4` row (no env var, prototype both times).
Expected: every process `CACHED` on the second run. This is the default-unchanged guarantee.

- [ ] **Step 6: Commit**

The pipeline lives under `specs/`, which is untracked for this feature. Commit only if the
spec directory is later tracked; otherwise record the runbook results in the spec document.

```bash
git status --short specs/260916-task-hash-spec/
```

---

### Task 11: The `global/v1` spec (applied on the seqeralabs #47 branch)

This task **does not apply to master** — `GlobalTaskHasher` and `FileIdentityStrategy` exist only on `seqeralabs/nextflow` branch `global-cache-provenance` (PR #47). Do it after Tasks 1–10 are merged there, or on a branch that merges both.

**Files (on the #47 branch):**
- Create: `plugins/nf-cloudcache-global/src/main/nextflow/cache/global/GlobalHashSpec.groovy`
- Create: `plugins/nf-cloudcache-global/src/main/nextflow/cache/global/GlobalTaskHashSpecFactory.groovy`
- Delete: `plugins/nf-cloudcache-global/src/main/nextflow/cache/global/GlobalTaskHasher.groovy`
- Modify: `plugins/nf-cloudcache-global/src/resources/META-INF/extensions.idx`
- Test: `plugins/nf-cloudcache-global/src/test/nextflow/cache/global/GlobalHashSpecTest.groovy`

**Interfaces:**
- Consumes: `TaskHashSpec`, `KeyBinding`, `HashKey`, `Contributors`, `TaskHashSpecFactory` (Tasks 2–8); `FileIdentityStrategy` from the plugin.
- Produces: `GlobalHashSpec.of(FileIdentityStrategy identity)` returning `TaskHashSpec` with id `global/v1`.

- [ ] **Step 1: Write the failing test**

```groovy
// plugins/nf-cloudcache-global/src/test/nextflow/cache/global/GlobalHashSpecTest.groovy
package nextflow.cache.global

import nextflow.cache.global.identity.FileIdentityStrategy
import nextflow.processor.hash.HashKey
import nextflow.processor.hash.StdSpecs
import spock.lang.Specification

class GlobalHashSpecTest extends Specification {

    def 'global/v1 is std/v4 without the run-specific keys'() {
        given:
        def spec = GlobalHashSpec.of(Mock(FileIdentityStrategy))

        expect:
        spec.id == 'global/v1'
        and: 'the hash is portable across runs and processes'
        !(HashKey.SESSION_ID in spec.keys())
        !(HashKey.PROCESS_NAME in spec.keys())
        and: 'every other key of std/v4 is retained, in the same order'
        spec.keys() == StdSpecs.STD_V4.keys() - [HashKey.SESSION_ID, HashKey.PROCESS_NAME]
    }

    def 'the identity strategy is part of the spec identity'() {
        given:
        def a = GlobalHashSpec.of(Mock(FileIdentityStrategy) { canonicalName() >> 'sample' })
        def b = GlobalHashSpec.of(Mock(FileIdentityStrategy) { canonicalName() >> 'provenance' })

        expect:
        a.fingerprint() != b.fingerprint()
    }
}
```

- [ ] **Step 2: Run the test and confirm it fails**

Run: `./gradlew :plugins:nf-cloudcache-global:test --tests "nextflow.cache.global.GlobalHashSpecTest"`
Expected: FAIL — `unable to resolve class GlobalHashSpec`.

- [ ] **Step 3: Implement the spec**

`FileIdentityStrategy` must expose `canonicalName()`; add it if absent, returning a stable
short name per implementation (`'sample'`, `'provenance'`). Swapping strategies changes
hashes, so it must move the fingerprint.

```groovy
// plugins/nf-cloudcache-global/src/main/nextflow/cache/global/GlobalHashSpec.groovy
package nextflow.cache.global

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.cache.global.identity.FileIdentityStrategy
import nextflow.file.FileHolder
import nextflow.util.ArrayBag
import org.apache.commons.collections4.Bag
import nextflow.processor.hash.Contributor
import nextflow.processor.hash.Contributors
import nextflow.processor.hash.HashContext
import nextflow.processor.hash.HashKey
import nextflow.processor.hash.KeyBinding
import nextflow.processor.hash.TaskHashSpec
import nextflow.util.EncodingRules

/**
 * The content-only task hash spec used by the global cache: std/v4 minus the two
 * run-specific keys, with file inputs identified by content rather than by path.
 */
@CompileStatic
class GlobalHashSpec {

    static TaskHashSpec of(FileIdentityStrategy identity) {
        return new TaskHashSpec('global/v1', [
            new KeyBinding(HashKey.TASK_SOURCE, Contributors.TASK_SOURCE),
            new KeyBinding(HashKey.CONTAINER, Contributors.CONTAINER),
            new KeyBinding(HashKey.INPUTS, inputsByIdentity(identity)),
            new KeyBinding(HashKey.EVAL_OUTPUTS, Contributors.EVAL_OUTPUTS_RAW_MAP),
            new KeyBinding(HashKey.SCRIPT_VARS, Contributors.SCRIPT_VARS),
            new KeyBinding(HashKey.BIN_ENTRIES, Contributors.BIN_ENTRIES),
            new KeyBinding(HashKey.MODULE_BUNDLE, Contributors.MODULE_BUNDLE),
            new KeyBinding(HashKey.ENV_MODULES, Contributors.ENV_MODULES),
            new KeyBinding(HashKey.CONDA, Contributors.CONDA),
            new KeyBinding(HashKey.SPACK, Contributors.SPACK_AND_ARCH),
            new KeyBinding(HashKey.STUB_MARKER, Contributors.STUB_MARKER)
        ], EncodingRules.RECORD_TYPES)
    }

    /**
     * As Contributors.INPUTS_RAW, but each value passes through the identity strategy,
     * recursing into nested collections exactly as GlobalTaskHasher.contributeInput did.
     */
    private static Contributor inputsByIdentity(FileIdentityStrategy identity) {
        return Contributors.of("inputs.identity:${identity.canonicalName()}") { HashContext ctx ->
            final out = new ArrayList<Object>()
            for( final entry : ctx.task.inputs ) {
                out.add(entry.key.name)
                out.add(contributeInput(identity, entry.value))
            }
            return out
        }
    }

    /** GlobalTaskHasher.contributeInput, with the strategy passed rather than held. */
    private static Object contributeInput(FileIdentityStrategy identity, Object value) {
        if( identity == null )
            return value
        // -- the two file shapes this strategy can identify
        if( value instanceof FileHolder ) {
            final id = identity.identify((FileHolder)value)
            return id != null ? id : value
        }
        if( value instanceof Path ) {
            final id = identity.identify((Path)value)
            return id != null ? id : value
        }
        // -- containers, in HashBuilder's own order
        if( value instanceof byte[] )
            return value                                    // a leaf there, not a container
        if( value instanceof Object[] )
            return contributeAll(identity, Arrays.asList((Object[]) value)).toArray()
        if( value instanceof Map )
            return contributeMap(identity, (Map) value)
        if( value instanceof Map.Entry ) {
            final e = (Map.Entry) value
            return Map.entry(contributeInput(identity, e.getKey()), contributeInput(identity, e.getValue()))
        }
        if( value instanceof Bag || value instanceof Set )
            return new ArrayBag<Object>(contributeAll(identity, (Collection) value))
        if( value instanceof Collection )
            return contributeAll(identity, (Collection) value)
        return value
    }

    private static List<Object> contributeAll(FileIdentityStrategy identity, Collection<?> items) {
        final list = items instanceof List ? (List) items : new ArrayList<Object>(items)
        final ids = identity.identifyAll((List<Object>) list)
        final out = new ArrayList<Object>(list.size())
        for( int i=0; i<list.size(); i++ ) {
            final item = list.get(i)
            final id = ids != null && i < ids.size() ? ids.get(i) : null
            if( id != null )
                out.add(id)
            else if( item instanceof FileHolder || item instanceof Path )
                out.add(item)                               // not identifiable -> default hashing
            else
                out.add(contributeInput(identity, item))    // a nested container, or a non-file value
        }
        return out
    }

    private static Map<Object,Object> contributeMap(FileIdentityStrategy identity, Map<?,?> value) {
        // Keys and values are flattened into ONE batch, so a map carrying files is identified on the
        // pool like any other collection instead of one entry at a time -- `contributeAll` hands the
        // whole list to the strategy, which may spread it.
        final flat = new ArrayList<Object>(value.size() * 2)
        for( Map.Entry<?,?> entry : value.entrySet() ) {
            flat.add(entry.getKey())
            flat.add(entry.getValue())
        }
        final done = contributeAll(identity, flat)
        final out = new LinkedHashMap<Object,Object>(value.size())
        for( int i=0; i<done.size(); i+=2 )
            out.put(done.get(i), done.get(i+1))
        return out.size() == value.size() ? out : (Map<Object,Object>) value
    }
}
```

```groovy
// plugins/nf-cloudcache-global/src/main/nextflow/cache/global/GlobalTaskHashSpecFactory.groovy
package nextflow.cache.global

import groovy.transform.CompileStatic
import nextflow.processor.TaskRun
import nextflow.processor.hash.TaskHashSpec
import nextflow.processor.hash.TaskHashSpecFactory
import org.pf4j.Extension

@Extension
@CompileStatic
class GlobalTaskHashSpecFactory implements TaskHashSpecFactory {

    @Override
    TaskHashSpec create(TaskRun task) {
        final config = GlobalCacheConfig.of(task.processor.session)
        if( !config.enabled ) {
            return null
        }
        return GlobalHashSpec.of(config.fileIdentityStrategy())
    }
}
```

- [ ] **Step 4: Retire GlobalTaskHasher**

The three `contribute*` methods above are `GlobalTaskHasher`'s, with the strategy passed as a
parameter rather than held as a field; diff them against the original to confirm nothing else
changed. Then delete `GlobalTaskHasher.groovy` and replace its entry in `extensions.idx` with
`GlobalTaskHashSpecFactory`. Verify the exact import for `Bag` against the original file — it is
whatever `GlobalTaskHasher` already imports.

- [ ] **Step 5: Run the plugin suite**

Run: `./gradlew :plugins:nf-cloudcache-global:test`
Expected: PASS, including the plugin's existing hashing tests.

- [ ] **Step 6: Run the global half of the runbook**

```bash
rm -rf work .nextflow*
<#47 build> run main.nf -profile withtooling,globalcache
<prototype+#47 build> run main.nf -profile withtooling,globalcache -resume
```
Expected: every process `CACHED`.

- [ ] **Step 7: Commit**

```bash
git add plugins/nf-cloudcache-global/src/main/nextflow/cache/global/GlobalHashSpec.groovy \
        plugins/nf-cloudcache-global/src/main/nextflow/cache/global/GlobalTaskHashSpecFactory.groovy \
        plugins/nf-cloudcache-global/src/test/nextflow/cache/global/GlobalHashSpecTest.groovy \
        plugins/nf-cloudcache-global/src/resources/META-INF/extensions.idx
git rm plugins/nf-cloudcache-global/src/main/nextflow/cache/global/GlobalTaskHasher.groovy
git commit -s -m "Express the global cache hasher as a task hash spec"
```
