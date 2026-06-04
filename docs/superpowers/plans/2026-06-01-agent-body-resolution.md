# Agent Body Resolution Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the full directive-rich `agent` DSL form (`model`, `instruction`, `tools`, `maxIterations`, typed `input:`/`output:`, and a `${var}`-interpolated `prompt:`) pass name resolution under typed mode and load cleanly via `ScriptLoaderV2`.

**Architecture:** Mirror the existing `process` resolution path. Introduce an `AgentDsl` DSL-scope interface (parallel to `ProcessDsl`) that declares the four agent directive method shapes, then override `VariableScopeVisitor.visitAgent` to (a) push the `AgentDsl` definition scope, (b) declare the agent's `input:` parameters as locals so the `prompt:` template can reference them, (c) push the `AgentDsl.DirectiveDsl` scope while visiting directives, and (d) visit the prompt and typed outputs. No runtime/execution changes — execution remains stubbed (`AgentDef.run()` throws). Tool-reference resolution (`tools fastqc`) stays out of scope; `tools()` is exercised empty.

**Tech Stack:** Java (nf-lang front-end / ANTLR-driven AST + Groovy AST visitors), Groovy + Spock tests, Gradle.

---

## Background / why this plan exists

The `agent` keyword already parses, lowers to Groovy (`agent('name', {})`), and registers an `AgentDef` at runtime. But the **body** does not resolve under typed mode:

- `VariableScopeVisitor` has **no** `visitAgent` override. It inherits `ScriptVisitorSupport.visitAgent`, which visits `directives`, `outputs`, and `prompt` but never pushes a DSL scope and never declares the `input:` parameters. As a result:
  - The `prompt:` template `${question}` references an undeclared variable → `DynamicVariablesVisitor` (in `ScriptResolveVisitor`) reports `` `question` is not defined``.
  - The directive method calls (`model`, `instruction`, `tools`, `maxIterations`) are checked against the default `ScriptDsl` scope, where they do not exist.

The existing `AgentParserTest` cases mask this because their `parse()` helper asserts only on `Phases.SYNTAX` errors (`TestUtils.hasSyntaxErrors` filters to `PhaseAware.getPhase() == Phases.SYNTAX`). `NAME_RESOLUTION` errors are silently ignored there. `TestUtils.check(...)` returns **all** errors regardless of phase — that is the harness this plan uses to drive the failing test.

### Key source references (read before starting)

- `modules/nf-lang/src/main/java/nextflow/script/control/VariableScopeVisitor.java`
  - `visitProcessV2(ProcessNodeV2)` at lines ~339–374 — the pattern to mirror.
  - `visitDirectives(Statement, String, boolean)` at ~446 — reused as-is.
  - `visitTypedOutputs(Statement, String)` at ~311 — reused as-is.
  - `asFlatParams(...)` helper (used at line 345) — reused as-is.
  - `currentDefinition` field (set/cleared around process visits) — reused.
  - There is currently **no** `visitAgent` here (only `declareMethod(agentNode)` at line 122 and the `AgentNode` branch in `dataflowMethodType` at ~764).
- `modules/nf-lang/src/main/java/nextflow/script/control/ScriptResolveVisitor.java`
  - `visitAgent(AgentNode)` at ~148–156 already resolves types for inputs/directives/outputs/prompt. **No change needed** — it depends on the scopes that `VariableScopeVisitor` builds.
- `modules/nf-lang/src/main/java/nextflow/script/dsl/ProcessDsl.java` — the DSL-scope template (interface extends `DslScope`, nested `DirectiveDsl`, `@Description` Java text blocks; `Description`/`Constant` live in the same package so no import is needed).
- `modules/nf-lang/src/main/java/nextflow/script/ast/AgentNode.java` — fields: `Statement directives`, `Parameter[] inputs`, `Statement outputs`, `Statement prompt`.
- `modules/nf-lang/src/main/java/nextflow/script/dsl/DslScope.java` — empty marker interface.

---

## File structure

### Files to create

- `modules/nf-lang/src/main/java/nextflow/script/dsl/AgentDsl.java` — DSL scope for agent definitions. Outer interface (definition scope where inputs are declared) plus a nested `DirectiveDsl` interface declaring the four agent directives. One responsibility: describe the agent body's resolvable method surface.

### Files to modify

- `modules/nf-lang/src/main/java/nextflow/script/control/VariableScopeVisitor.java` — add an `import nextflow.script.dsl.AgentDsl;` and an `@Override public void visitAgent(AgentNode node)` method that builds the agent body scopes.

### Test files

- `modules/nf-lang/src/test/groovy/nextflow/script/parser/AgentParserTest.groovy` — add a resolution test using `check(...)` asserting **no** errors for the directive-rich agent (drives the front-end fix).
- `modules/nextflow/src/test/groovy/nextflow/script/parser/v2/AgentScriptLoadingTest.groovy` — add an end-to-end test that loads a directive-rich agent via `ScriptLoaderV2` and asserts the `AgentDef` registers (proves the full load path works, not just `nf-lang` analysis).

---

## Pre-flight

- [ ] **Step 0.1: Confirm baseline tests pass**

Run:
```bash
./gradlew :nf-lang:test --tests "nextflow.script.parser.AgentParserTest"
```
Expected: PASS (3 tests). These are the pre-existing parser tests; they must stay green.

- [ ] **Step 0.2: Confirm the agent runtime/loading tests pass**

Run:
```bash
./gradlew :nextflow:test --tests "nextflow.script.parser.v2.AgentScriptLoadingTest"
```
Expected: PASS (1 test — the minimal `agent foo { prompt: "..." }` fixture).

---

## Task 1: Failing resolution test for the directive-rich agent body

**Files:**
- Test: `modules/nf-lang/src/test/groovy/nextflow/script/parser/AgentParserTest.groovy`

- [ ] **Step 1.1: Add the failing test**

Append this method to `AgentParserTest` (before the closing brace). It uses `check(...)`, which returns **all** errors (not just syntax), so the current `NAME_RESOLUTION` gap surfaces:

```groovy
    def 'should resolve directives and prompt variables in an agent body'() {
        when:
        def errors = check('''\
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                    question: String

                output:
                    plan: String

                prompt:
                """
                Question: ${question}
                """
            }
            ''')

        then:
        errors.isEmpty()
    }
```

- [ ] **Step 1.2: Run the test to verify it fails**

Run:
```bash
./gradlew :nf-lang:test --tests "nextflow.script.parser.AgentParserTest"
```
Expected: FAIL on `should resolve directives and prompt variables in an agent body`. The failure lists one or more `NAME_RESOLUTION` errors — at minimum `` `question` is not defined`` (the prompt variable), and unresolved agent directive method(s). The other three `AgentParserTest` cases still PASS.

> If this test unexpectedly PASSES, stop: the assumption behind this plan (that the agent body does not resolve) is wrong. Re-read `VariableScopeVisitor` for a `visitAgent` override and reconcile before continuing.

- [ ] **Step 1.3: Commit the failing test**

```bash
git add modules/nf-lang/src/test/groovy/nextflow/script/parser/AgentParserTest.groovy
git commit -s -m "test: failing resolution test for directive-rich agent body"
```

---

## Task 2: Create the `AgentDsl` DSL scope

**Files:**
- Create: `modules/nf-lang/src/main/java/nextflow/script/dsl/AgentDsl.java`

- [ ] **Step 2.1: Write `AgentDsl.java`**

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
package nextflow.script.dsl;

/**
 * DSL scope for agent definitions.
 *
 * Mirrors {@link ProcessDsl}: the outer interface is the definition scope
 * (where the agent's typed `input:` parameters are declared as locals), and
 * the nested {@link DirectiveDsl} declares the agent directive methods that
 * may appear at the top of an agent body.
 */
public interface AgentDsl extends DslScope {

    interface DirectiveDsl extends DslScope {

        @Description("""
            The `model` directive selects the LLM in `provider/model` form (e.g. `openai/gpt-5-mini`).
        """)
        void model(String value);

        @Description("""
            The `instruction` directive sets the agent system prompt (its role/persona).
        """)
        void instruction(String value);

        @Description("""
            The `tools` directive declares the modules the agent may invoke as tools.
        """)
        void tools(Object... values);

        @Description("""
            The `maxIterations` directive caps the LLM tool-calling loop.
        """)
        void maxIterations(Integer value);
    }
}
```

- [ ] **Step 2.2: Verify it compiles**

Run:
```bash
./gradlew :nf-lang:compileJava
```
Expected: BUILD SUCCESSFUL.

- [ ] **Step 2.3: Commit**

```bash
git add modules/nf-lang/src/main/java/nextflow/script/dsl/AgentDsl.java
git commit -s -m "feat: add AgentDsl directive scope"
```

---

## Task 3: Override `VariableScopeVisitor.visitAgent`

**Files:**
- Modify: `modules/nf-lang/src/main/java/nextflow/script/control/VariableScopeVisitor.java`

- [ ] **Step 3.1: Add the `AgentDsl` import**

Find the existing DSL imports block (it contains `import nextflow.script.dsl.ProcessDsl;` around line 45). Add immediately above it:

```java
import nextflow.script.dsl.AgentDsl;
```

(Keep imports alphabetically ordered within the `nextflow.script.dsl` group: `AgentDsl` precedes `Constant`.)

- [ ] **Step 3.2: Add the `visitAgent` override**

Insert this method immediately **before** `visitProcessV2` (currently at line ~339). It mirrors `visitProcessV2`: push the definition scope, declare inputs as locals (so the prompt can reference them), visit directives under the directive scope, then visit prompt and typed outputs.

```java
    @Override
    public void visitAgent(AgentNode node) {
        vsc.pushScope(AgentDsl.class);
        currentDefinition = node;
        node.setVariableScope(currentScope());

        for( var input : asFlatParams(node.inputs) ) {
            vsc.declare(input, input);

            // suppress "unused variable" warnings for inputs only referenced in the prompt
            vsc.findVariableDeclaration(input.getName(), input);
        }

        vsc.pushScope(AgentDsl.DirectiveDsl.class);
        visitDirectives(node.directives, "agent directive", false);
        vsc.popScope();

        // the prompt template may reference input parameters
        visit(node.prompt);

        visitTypedOutputs(node.outputs, "Agent output");

        currentDefinition = null;
        vsc.popScope();
    }
```

- [ ] **Step 3.3: Run the resolution test to verify it passes**

Run:
```bash
./gradlew :nf-lang:test --tests "nextflow.script.parser.AgentParserTest"
```
Expected: PASS (4 tests — including `should resolve directives and prompt variables in an agent body`).

- [ ] **Step 3.4: Commit**

```bash
git add modules/nf-lang/src/main/java/nextflow/script/control/VariableScopeVisitor.java
git commit -s -m "feat: resolve agent body directives, inputs and prompt vars"
```

---

## Task 4: End-to-end load test for a directive-rich agent

**Files:**
- Test: `modules/nextflow/src/test/groovy/nextflow/script/parser/v2/AgentScriptLoadingTest.groovy`

- [ ] **Step 4.1: Add the directive-rich load test**

Append this method to `AgentScriptLoadingTest` (before the closing brace). It exercises the full `ScriptLoaderV2` path — name resolution, Groovy lowering, and `runScript()` registration — for the directive-rich form. Execution is still stubbed (the lowered closure is empty per the current `ScriptToGroovyVisitor.visitAgent`), so registration succeeds without invoking the agent:

```groovy
    def 'should load a script with a directive-rich agent definition'() {
        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = '''
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                    question: String

                output:
                    plan: String

                prompt:
                """
                Question: ${question}
                """
            }

            workflow {
            }
            '''.stripIndent()

        when:
        parser.parse(file)
        parser.runScript()

        then:
        def definitions = ScriptMeta.get(parser.script).getDefinitions()
        definitions.any { it instanceof AgentDef && it.name == 'eval_agent' }

        cleanup:
        file.parent.deleteDir()
    }
```

- [ ] **Step 4.2: Run the loading tests**

Run:
```bash
./gradlew :nextflow:test --tests "nextflow.script.parser.v2.AgentScriptLoadingTest"
```
Expected: PASS (2 tests — the original minimal fixture plus the new directive-rich one).

> If `parse(file)` throws a compilation error mentioning an unresolved variable or method, the front-end fix from Task 3 did not take effect for the `:nextflow` module's view of `nf-lang`. Re-run `./gradlew :nextflow:compileTestGroovy` to force a rebuild of the `nf-lang` dependency, then re-run.

- [ ] **Step 4.3: Commit**

```bash
git add modules/nextflow/src/test/groovy/nextflow/script/parser/v2/AgentScriptLoadingTest.groovy
git commit -s -m "test: end-to-end load of a directive-rich agent definition"
```

---

## Task 5: Full verification

- [ ] **Step 5.1: Run the full nf-lang and agent suites**

Run:
```bash
./gradlew :nf-lang:test :nextflow:test --tests "*Agent*"
```
Expected: PASS. Specifically `AgentParserTest` (4), `AgentScriptLoadingTest` (2), `AgentDefTest`, `AgentFactoryTest` all green, with no regressions in the broader `nf-lang` suite.

- [ ] **Step 5.2: Confirm no unrelated regressions in nf-lang resolution**

Run:
```bash
./gradlew :nf-lang:test --tests "nextflow.script.control.*"
```
Expected: PASS (the `VariableScopeVisitor` change must not affect process/workflow/function resolution).

- [ ] **Step 5.3: Update the ADR implementation status**

Edit `adr/20260505-llm-agent-primitive.md`: in "Known limitations / what's missing", mark item 1 (Agent body resolution) as resolved, and in "Next steps" mark "Plan B" delivered, pointing to this plan file. Keep items 2–6 (closure population, execution, tool bridge, config scope) as still-open.

- [ ] **Step 5.4: Commit the ADR update**

```bash
git add adr/20260505-llm-agent-primitive.md
git commit -s -m "docs: mark agent body resolution (Plan B) delivered"
```

---

## What this plan does NOT deliver

These remain future work (unchanged from the ADR's Next steps):

1. **Closure population (Plan C)** — `ScriptToGroovyVisitor.visitAgent` still emits `agent('name', {})` with an empty body. Directives/inputs/outputs/prompt live on the `AgentNode` but are not propagated into the runtime closure.
2. **Agent execution (Plan D)** — `AgentDef.run()` still throws `UnsupportedOperationException`.
3. **`nf-agent` plugin** — langchain4j runtime, chat-model factory, tool dispatcher.
4. **Tool-reference resolution** — `tools fastqc, multiqc` (module references as directive args) is intentionally out of scope; this plan only exercises empty `tools()`. Resolving module references as tool args is part of the tool-bridge work.
5. **`agent` config scope** — `nextflow.config` `agent { ... }` parsing.
6. **Inheriting standard process directives** — `errorStrategy`, `maxRetries`, `time`, `cache` on agents (spec §2) are not added to `AgentDsl.DirectiveDsl` here.

---

## Self-Review

- **Spec coverage:** This plan covers exactly the spec §2 "Directives" table for the agent-specific four (`model`, `instruction`, `tools`, `maxIterations`) plus typed `input:`/`output:` and `${var}` prompt interpolation, at the resolution/load level. Execution (§3–§6) and the `nf-agent` plugin (§7) are explicitly out of scope and listed above. Tool-reference resolution and inherited process directives are flagged as deferred.
- **Placeholder scan:** No TBD/TODO/"handle edge cases" — every code step contains complete content.
- **Type consistency:** `AgentDsl.DirectiveDsl` is referenced with that exact nested name in `VariableScopeVisitor.visitAgent` (Step 3.2) and defined with that exact name in Step 2.1. `maxIterations(Integer)` matches the `maxIterations 20` literal. `visitTypedOutputs`, `visitDirectives`, `asFlatParams`, `currentDefinition`, `vsc.pushScope(Class)`, `vsc.declare(Variable, ASTNode)`, `vsc.findVariableDeclaration(String, ASTNode)` all match existing signatures in `VariableScopeVisitor`/`VariableScopeChecker`.
