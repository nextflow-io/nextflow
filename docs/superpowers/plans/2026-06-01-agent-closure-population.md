# Agent Closure Population Implementation Plan (Plan C)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the lowered `agent` Groovy closure carry the agent's directives, typed inputs, typed outputs, and prompt, and have the runtime build an `AgentDef` populated with that content (instead of the current empty-closure stub). Execution stays stubbed — `AgentDef.run()` continues to throw — but the populated structure is what Plan D's runner will consume.

**Architecture:** Mirror the `process` lowering+runtime pattern, collapsed to the simpler agent shape. A new `AgentToGroovyVisitorV2` lowers `AgentNode` into `agent('name', { <directives>; _input_(...); _output_(...); new PromptDef(...) })`. At runtime `BaseScript.agent(name, closure)` runs the closure against an `AgentBuilder` delegate (directives captured via `methodMissing` against a whitelist; inputs/outputs via `_input_`/`_output_`), takes the returned `PromptDef`, and builds an `AgentDef` exposing `model`/`instruction`/`tools`/`maxIterations`/`inputs`/`outputs`/`prompt`.

**Tech Stack:** Java (nf-lang AST→Groovy lowering), Groovy (runtime model), Spock, Gradle.

---

## Background (read first)

The current state after Plan B:
- `ScriptToGroovyVisitor.visitAgent` (`modules/nf-lang/src/main/java/nextflow/script/control/ScriptToGroovyVisitor.java:288-294`) emits `agent('name', {})` — an EMPTY closure. The directives/inputs/outputs/prompt live on the `AgentNode` but never reach runtime.
- `BaseScript.agent(String name, Closure body)` (`modules/nextflow/src/main/groovy/nextflow/script/BaseScript.groovy:140-144`) creates an `AgentFactory` and registers a bare `AgentDef` holding the (empty) closure.
- `AgentDef` (`modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy`) stores only `owner`/`name`/`simpleName`/`body`; `run()` throws `UnsupportedOperationException`.
- `AgentFactory` (`modules/nextflow/src/main/groovy/nextflow/script/AgentFactory.groovy`) is a thin `newAgent(name, body)` wrapper.

The pattern to mirror (process path, for reference — do NOT modify these):
- Lowering: `ProcessToGroovyVisitorV2.transform()` builds a multi-statement closure ending in a `new BodyDef(...)` expression; `ScriptToGroovyVisitor.visitProcessV2` delegates to it.
- Runtime: `BaseScript.processV2(name, closure)` clones the closure, sets a `ProcessDslV2` delegate with `DELEGATE_FIRST`, calls it to get the `BodyDef`, then `dsl.withBody(body).build()`.
- Directive capture: `ProcessBuilder.methodMissing(name, args)` validates the name and puts it into a `ProcessConfig` map.
- Body capture: `BodyDef(Closure closure, String source, String section)` holds the script closure + source text.

Agents are simpler than processes: inputs are typed `val` declarations only (no staging/env/files), there is exactly one output kind (a typed `val`), and the "body" is a `prompt:` template rather than a `script:`. So we collapse the 3-layer process design (Builder+Config+DslV2) into a single `AgentBuilder` delegate.

## Naming note

`nextflow.script.dsl.AgentDsl` already exists (the **compile-time** resolution scope interface from Plan B). The **runtime** delegate introduced here is a different type named `AgentBuilder` in package `nextflow.script`. Do not confuse or merge them.

---

## File structure

### Files to create
- `modules/nextflow/src/main/groovy/nextflow/script/PromptDef.groovy` — captures the prompt closure + source text (parallel to `BodyDef`, minimal).
- `modules/nextflow/src/main/groovy/nextflow/script/AgentBuilder.groovy` — runtime delegate: captures directives (`methodMissing` against a whitelist), `_input_`/`_output_`, holds inputs/outputs, `build()` → `AgentDef`. Includes small static-typed holders `AgentInput`/`AgentOutput` as nested classes.
- `modules/nf-lang/src/main/java/nextflow/script/control/AgentToGroovyVisitorV2.java` — lowers `AgentNode` to the populated `agent(...)` call.

### Files to modify
- `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` — store `config`/`inputs`/`outputs`/`prompt`; expose typed getters. `run()` still throws.
- `modules/nextflow/src/main/groovy/nextflow/script/BaseScript.groovy` — `agent(name, closure)` runs the closure against `AgentBuilder` and builds the populated `AgentDef`.
- `modules/nf-lang/src/main/java/nextflow/script/control/ScriptToGroovyVisitor.java` — `visitAgent` delegates to `AgentToGroovyVisitorV2`.
- `modules/nextflow/src/main/groovy/nextflow/script/AgentFactory.groovy` — DELETE (superseded by `AgentBuilder`) **only if** no other references remain; otherwise leave untouched and note it. Check with grep first.

### Test files
- `modules/nextflow/src/test/groovy/nextflow/script/AgentBuilderTest.groovy` — unit tests for directive capture, input/output capture, prompt, and `build()`.
- `modules/nextflow/src/test/groovy/nextflow/script/parser/v2/AgentScriptLoadingTest.groovy` — extend the directive-rich load test to assert the populated `AgentDef` content.

---

## Task 1: Runtime model — `PromptDef`, `AgentBuilder`, populated `AgentDef`

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/script/PromptDef.groovy`
- Create: `modules/nextflow/src/main/groovy/nextflow/script/AgentBuilder.groovy`
- Modify: `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/script/AgentBuilderTest.groovy`

- [ ] **Step 1.1: Write the failing unit test**

Create `modules/nextflow/src/test/groovy/nextflow/script/AgentBuilderTest.groovy`:

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.script

import spock.lang.Specification

class AgentBuilderTest extends Specification {

    def 'should capture directives, inputs, outputs and prompt then build an AgentDef'() {
        given:
        def builder = new AgentBuilder(null, 'eval_agent')

        when:
        builder.model('openai/gpt-5-mini')
        builder.instruction('You are helpful.')
        builder.tools()
        builder.maxIterations(20)
        builder._input_('question', String)
        builder._output_('plan', String)
        def prompt = new PromptDef({ "Question: x" }, 'Question: ${question}')
        def agent = builder.withPrompt(prompt).build()

        then:
        agent instanceof AgentDef
        agent.name == 'eval_agent'
        agent.model == 'openai/gpt-5-mini'
        agent.instruction == 'You are helpful.'
        agent.maxIterations == 20
        agent.tools == []
        agent.inputs*.name == ['question']
        agent.outputs*.name == ['plan']
        agent.prompt.source == 'Question: ${question}'
    }

    def 'should reject an unknown directive'() {
        given:
        def builder = new AgentBuilder(null, 'a')

        when:
        builder.bogusDirective('x')

        then:
        thrown(Exception)
    }

    def 'should fail to build without a prompt'() {
        given:
        def builder = new AgentBuilder(null, 'a')

        when:
        builder.build()

        then:
        thrown(IllegalStateException)
    }
}
```

- [ ] **Step 1.2: Run it to verify it fails to compile**

Run:
```bash
./gradlew :nextflow:test --tests "nextflow.script.AgentBuilderTest"
```
Expected: FAIL — compilation error (`AgentBuilder`/`PromptDef` do not exist yet).

- [ ] **Step 1.3: Write `PromptDef.groovy`**

Create `modules/nextflow/src/main/groovy/nextflow/script/PromptDef.groovy`:

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.script

import groovy.transform.CompileStatic

/**
 * Models the `prompt:` block of an agent definition. Mirrors {@link BodyDef}
 * but minimal: the prompt template is captured as a closure (evaluated per
 * invocation with the agent inputs in scope) plus its source text.
 */
@CompileStatic
class PromptDef implements Cloneable {

    final Closure closure
    final String source

    PromptDef(Closure closure, String source) {
        this.closure = closure
        this.source = source
    }

    @Override
    PromptDef clone() {
        (PromptDef) super.clone()
    }
}
```

- [ ] **Step 1.4: Write `AgentBuilder.groovy`**

Create `modules/nextflow/src/main/groovy/nextflow/script/AgentBuilder.groovy`. Model `methodMissing` on `ProcessBuilder` (verify its exact shape in `ProcessBuilder.groovy` before finalizing):

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.script

import groovy.transform.CompileStatic
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j

/**
 * Runtime delegate that captures an agent body when the lowered agent closure
 * runs. Directives are captured via {@link #methodMissing} against a known set;
 * inputs/outputs via {@code _input_}/{@code _output_}; the prompt via
 * {@link #withPrompt}. {@link #build} produces the populated {@link AgentDef}.
 *
 * Distinct from {@code nextflow.script.dsl.AgentDsl} (the compile-time
 * resolution scope).
 */
@Slf4j
@CompileStatic
class AgentBuilder {

    static final List<String> DIRECTIVES = ['model', 'instruction', 'tools', 'maxIterations']

    private BaseScript ownerScript
    private String agentName

    private final Map<String,Object> directives = new LinkedHashMap<>()
    private final List<AgentInput> inputs = new ArrayList<>()
    private final List<AgentOutput> outputs = new ArrayList<>()
    private PromptDef prompt

    AgentBuilder(BaseScript ownerScript, String agentName) {
        this.ownerScript = ownerScript
        this.agentName = agentName
    }

    protected void checkName(String name) {
        if( !DIRECTIVES.contains(name) )
            throw new IllegalArgumentException("Unknown agent directive `${name}`")
    }

    @PackageScope
    Object methodMissing(String name, Object args) {
        checkName(name)
        final values = args instanceof Object[] ? (args as List) : [args]
        directives.put(name, values.size() == 1 ? values[0] : values)
        return null
    }

    void _input_(String name, Class type) {
        inputs.add(new AgentInput(name, type))
    }

    void _output_(String name, Class type) {
        outputs.add(new AgentOutput(name, type))
    }

    AgentBuilder withPrompt(PromptDef prompt) {
        this.prompt = prompt
        return this
    }

    AgentDef build() {
        if( prompt == null )
            throw new IllegalStateException("Missing prompt in agent `${agentName}` definition")
        return new AgentDef(ownerScript, agentName, directives, inputs, outputs, prompt)
    }

    @CompileStatic
    static class AgentInput {
        final String name
        final Class type
        AgentInput(String name, Class type) { this.name = name; this.type = type }
    }

    @CompileStatic
    static class AgentOutput {
        final String name
        final Class type
        AgentOutput(String name, Class type) { this.name = name; this.type = type }
    }
}
```

- [ ] **Step 1.5: Rewrite `AgentDef.groovy` to hold the populated structure**

Replace the body of `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` with (keep the license header that's already there; keep `extends BindableDef implements ChainableDef` and the existing `cloneWithName` shape, adapting fields):

```groovy
package nextflow.script

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.script.AgentBuilder.AgentInput
import nextflow.script.AgentBuilder.AgentOutput

/**
 * Runtime model for an agent definition. Holds the captured directives,
 * inputs, outputs and prompt. Execution is delegated to the agent runner
 * (nf-agent plugin) and is not yet implemented — {@link #run} throws.
 */
@Slf4j
@CompileStatic
class AgentDef extends BindableDef implements ChainableDef {

    static final String TYPE = 'agent'

    private BaseScript owner
    private String name
    private String simpleName
    private Map<String,Object> directives
    private List<AgentInput> inputs
    private List<AgentOutput> outputs
    private PromptDef prompt

    AgentDef(BaseScript owner, String name, Map<String,Object> directives, List<AgentInput> inputs, List<AgentOutput> outputs, PromptDef prompt) {
        this.owner = owner
        this.name = name
        this.simpleName = name
        this.directives = directives
        this.inputs = inputs
        this.outputs = outputs
        this.prompt = prompt
    }

    @Override String getType() { TYPE }
    @Override String getName() { name }
    String getSimpleName() { simpleName }
    BaseScript getOwner() { owner }

    String getModel() { directives.get('model') as String }
    String getInstruction() { directives.get('instruction') as String }
    List getTools() { (directives.get('tools') ?: []) as List }
    Integer getMaxIterations() { directives.get('maxIterations') as Integer }
    List<AgentInput> getInputs() { inputs }
    List<AgentOutput> getOutputs() { outputs }
    PromptDef getPrompt() { prompt }

    @Override
    ComponentDef cloneWithName(String name) {
        def copy = (AgentDef) this.clone()
        copy.@name = name
        copy.@simpleName = name.contains(':') ? name.tokenize(':').last() : name
        return copy
    }

    @Override
    Object run(Object[] args) {
        throw new UnsupportedOperationException('agent execution not yet implemented')
    }
}
```

> Before finalizing: check whether `BindableDef`/`ChainableDef` require any other abstract methods (e.g. `getDeclaredInputs`/arity). The previous `AgentDef` compiled with just `getType`/`getName`/`cloneWithName`/`run`, so this superset should too. If the compiler demands more, mirror what `ProcessDef` returns.

- [ ] **Step 1.6: Run the unit test to verify it passes**

Run:
```bash
./gradlew :nextflow:test --tests "nextflow.script.AgentBuilderTest"
```
Expected: PASS (3 tests).

- [ ] **Step 1.7: Commit**

```bash
git add modules/nextflow/src/main/groovy/nextflow/script/PromptDef.groovy \
        modules/nextflow/src/main/groovy/nextflow/script/AgentBuilder.groovy \
        modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy \
        modules/nextflow/src/test/groovy/nextflow/script/AgentBuilderTest.groovy
git commit -s -m "feat: populate AgentDef via AgentBuilder runtime delegate"
```

---

## Task 2: Lowering — `AgentToGroovyVisitorV2` + wire `BaseScript.agent`

**Files:**
- Create: `modules/nf-lang/src/main/java/nextflow/script/control/AgentToGroovyVisitorV2.java`
- Modify: `modules/nf-lang/src/main/java/nextflow/script/control/ScriptToGroovyVisitor.java`
- Modify: `modules/nextflow/src/main/groovy/nextflow/script/BaseScript.groovy`
- Modify (delete if unreferenced): `modules/nextflow/src/main/groovy/nextflow/script/AgentFactory.groovy`
- Test: `modules/nextflow/src/test/groovy/nextflow/script/parser/v2/AgentScriptLoadingTest.groovy`

- [ ] **Step 2.1: Extend the end-to-end load test to assert populated content**

In `AgentScriptLoadingTest.groovy`, change the `then:` of `'should load a script with a directive-rich agent definition'` to inspect the `AgentDef`:

```groovy
        then:
        def definitions = ScriptMeta.get(parser.script).getDefinitions()
        def agent = definitions.find { it instanceof AgentDef && it.name == 'eval_agent' } as AgentDef
        agent != null
        agent.model == 'openai/gpt-5-mini'
        agent.instruction == 'You are helpful.'
        agent.maxIterations == 20
        agent.tools == []
        agent.inputs*.name == ['question']
        agent.outputs*.name == ['plan']
        agent.prompt != null
        agent.prompt.source.contains('Question:')
```

- [ ] **Step 2.2: Run it to verify it fails**

Run:
```bash
./gradlew :nextflow:test --tests "nextflow.script.parser.v2.AgentScriptLoadingTest"
```
Expected: FAIL — the agent loads but its `model`/`inputs`/etc. are empty/null because the lowered closure is still empty.

- [ ] **Step 2.3: Write `AgentToGroovyVisitorV2.java`**

Create `modules/nf-lang/src/main/java/nextflow/script/control/AgentToGroovyVisitorV2.java`. Use `ProcessToGroovyVisitorV2` as the structural reference for imports and helpers (`closureX`, `block`, `stmt`, `callThisX`, `constX`, `classX`, `createX`, `args`, `ScriptToGroovyHelper`). The emitted closure must end in a `new PromptDef(...)` expression so the runtime delegate returns it.

```java
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package nextflow.script.control;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import nextflow.script.ast.AgentNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.VariableScope;
import org.codehaus.groovy.ast.expr.Expression;
import org.codehaus.groovy.ast.stmt.ExpressionStatement;
import org.codehaus.groovy.ast.stmt.Statement;
import org.codehaus.groovy.control.SourceUnit;

import static nextflow.script.ast.ASTUtils.*;
import static org.codehaus.groovy.ast.tools.GeneralUtils.*;

/**
 * Lowers an {@link AgentNode} to a runtime {@code agent('name', { ... })} call.
 * The generated closure carries the directives, typed inputs/outputs and a
 * {@code PromptDef}, mirroring how {@code ProcessToGroovyVisitorV2} lowers a
 * process body.
 */
public class AgentToGroovyVisitorV2 {

    private SourceUnit sourceUnit;

    private ScriptToGroovyHelper sgh;

    public AgentToGroovyVisitorV2(SourceUnit sourceUnit) {
        this.sourceUnit = sourceUnit;
        this.sgh = new ScriptToGroovyHelper(sourceUnit);
    }

    public Statement transform(AgentNode node) {
        var statements = new ArrayList<Statement>();
        statements.add(node.directives);
        statements.add(agentInputs(node.inputs));
        statements.add(agentOutputs(node.outputs));
        statements.add(agentPrompt(node.prompt));
        var body = closureX(block(new VariableScope(), statements));
        return stmt(callThisX("agent", args(constX(node.getName()), body)));
    }

    private Statement agentInputs(Parameter[] inputs) {
        var statements = Arrays.stream(inputs)
            .map(input -> stmt(callThisX("_input_", args(constX(input.getName()), classX(input.getType())))))
            .map(s -> (Statement) s)
            .toList();
        return block(null, statements);
    }

    private Statement agentOutputs(Statement outputs) {
        var statements = new ArrayList<Statement>();
        for( var stmt : asBlockStatements(outputs) ) {
            var output = ((ExpressionStatement) stmt).getExpression();
            var target = targetOf(output);
            if( target != null )
                statements.add(stmt(callThisX("_output_", args(constX(target.getName()), classX(target.getType())))));
        }
        return block(null, statements);
    }

    private org.codehaus.groovy.ast.expr.VariableExpression targetOf(Expression output) {
        if( output instanceof org.codehaus.groovy.ast.expr.VariableExpression ve )
            return ve;
        if( output instanceof nextflow.script.ast.AssignmentExpression ae && ae.getLeftExpression() instanceof org.codehaus.groovy.ast.expr.VariableExpression ve )
            return ve;
        return null;
    }

    private Statement agentPrompt(Statement prompt) {
        var expr = ((ExpressionStatement) prompt).getExpression();
        return stmt(createX(
            "nextflow.script.PromptDef",
            args(
                closureX(stmt(expr)),
                constX(sgh.getSourceText(prompt))
            )
        ));
    }
}
```

> Verify the exact helper names/imports against `ProcessToGroovyVisitorV2.java`: `block`, `closureX`, `stmt`, `callThisX`, `constX`, `classX`, `createX`, `args`, `asBlockStatements`, and `ScriptToGroovyHelper.getSourceText`. Adjust import statics to match (the process visitor is the source of truth). The `targetOf` helper mirrors the output-target extraction in `ScriptResolveVisitor.resolveTypedOutputs`/`VariableScopeVisitor.visitTypedOutputs`. Confirm `AssignmentExpression`'s package (`nextflow.script.ast` per its use in `ScriptResolveVisitor`).

- [ ] **Step 2.4: Wire `ScriptToGroovyVisitor.visitAgent` to delegate**

Replace the body of `visitAgent` in `ScriptToGroovyVisitor.java` (currently lines ~288-294) with:

```java
    @Override
    public void visitAgent(AgentNode node) {
        checkReservedMethodName(node, "agent");
        var result = new AgentToGroovyVisitorV2(sourceUnit).transform(node);
        moduleNode.addStatement(result);
    }
```

(Mirror exactly how `visitProcessV2` delegates to `ProcessToGroovyVisitorV2`.)

- [ ] **Step 2.5: Update `BaseScript.agent` to run the closure against `AgentBuilder`**

Replace `BaseScript.agent` (`modules/nextflow/src/main/groovy/nextflow/script/BaseScript.groovy:140-144`) with a delegate-first execution mirroring `processV2`:

```groovy
    protected void agent(String name, Closure<PromptDef> body) {
        final builder = new AgentBuilder(this, name)
        final cl = (Closure<PromptDef>) body.clone()
        cl.setDelegate(builder)
        cl.setResolveStrategy(Closure.DELEGATE_FIRST)
        final prompt = cl.call()
        final agent = builder.withPrompt(prompt).build()
        meta.addDefinition(agent)
    }
```

Ensure `nextflow.script.PromptDef` and `nextflow.script.AgentBuilder` are resolvable (same package — no import needed).

- [ ] **Step 2.6: Remove the now-unused `AgentFactory` if nothing references it**

Run:
```bash
grep -rn "AgentFactory" modules/ plugins/ --include=*.groovy --include=*.java
```
If the only matches are the file itself and its test, delete both:
```bash
git rm modules/nextflow/src/main/groovy/nextflow/script/AgentFactory.groovy \
       modules/nextflow/src/test/groovy/nextflow/script/AgentFactoryTest.groovy
```
If other code references it, leave it and note this as a concern instead.

- [ ] **Step 2.7: Run the load test to verify it passes**

Run:
```bash
./gradlew :nextflow:test --tests "nextflow.script.parser.v2.AgentScriptLoadingTest"
```
Expected: PASS (2 tests) — the directive-rich agent now exposes `model`/`instruction`/`maxIterations`/`inputs`/`outputs`/`prompt`.

- [ ] **Step 2.8: Commit**

```bash
git add -A
git commit -s -m "feat: lower agent body into populated runtime closure"
```

---

## Task 3: Verification

- [ ] **Step 3.1: Run the full agent + lang suites**

Run:
```bash
./gradlew :nf-lang:test :nextflow:test --tests "*Agent*"
```
Expected: PASS — `AgentParserTest` (4), `AgentBuilderTest` (3), `AgentScriptLoadingTest` (2), `AgentDefTest` (adapt if it referenced the old constructor — see 3.2). No regressions.

- [ ] **Step 3.2: Reconcile `AgentDefTest`**

`AgentDefTest` may have constructed the old `AgentDef(owner, body, name)`. Update it to the new constructor (or to build via `AgentBuilder`) so it compiles and meaningfully tests the populated model. Run:
```bash
./gradlew :nextflow:test --tests "nextflow.script.AgentDefTest"
```
Expected: PASS.

- [ ] **Step 3.3: Regression guard on lowering**

Run:
```bash
./gradlew :nf-lang:test --tests "nextflow.script.control.*"
```
Expected: PASS (the `ScriptToGroovyVisitor` change must not affect process/workflow lowering).

- [ ] **Step 3.4: Commit any test reconciliation**

```bash
git add -A && git commit -s -m "test: reconcile agent tests with populated AgentDef"
```

---

## What this plan does NOT deliver
- **Execution (Plan D):** `AgentDef.run()` still throws. No langchain4j, no LLM calls, no tool dispatch.
- **Tool-reference resolution:** `tools fastqc` still unsupported; only empty `tools()` exercised. The lowering passes whatever args appear, but resolving module references is Plan D.
- **`agent` config scope, valRefs/caching on the prompt, structured outputs.**

## Self-Review
- **Spec coverage:** Covers ADR limitation #3 (empty lowered closure) end-to-end: directives, inputs, outputs, and prompt now reach an `AgentDef`. Execution stays out of scope (limitation #4 → Plan D).
- **Placeholder scan:** Complete code in every step; the two ">" notes are verification reminders, not deferred work.
- **Type consistency:** `AgentBuilder(BaseScript, String)`, `_input_(String, Class)`, `_output_(String, Class)`, `withPrompt(PromptDef)`, `build()→AgentDef`, `AgentDef(owner, name, Map, List<AgentInput>, List<AgentOutput>, PromptDef)`, `PromptDef(Closure, String)` — used identically across the lowering (`agent`, `_input_`, `_output_`, `new PromptDef(...)`), the runtime (`BaseScript.agent`), and the tests.
