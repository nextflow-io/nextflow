# Agent Syntax & Runtime Model — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the `agent` keyword to the Nextflow DSL and the runtime `AgentDef` / `AgentFactory` model, so agent definitions can be parsed, resolved, and loaded into a pipeline. Execution is intentionally out of scope for this plan — a future plan adds the `nf-agent` plugin and the langchain4j-driven runner.

**Architecture:** Mirror the existing `process` definition pipeline end-to-end. Lexer adds an `AGENT` keyword token, parser adds `agentDef`/`agentBody` rules, AST gains an `AgentNode`, the `ScriptVisitor` and `ScriptNode` collections are extended to carry agents, and the `ScriptAstBuilder` builds nodes from the parse tree. On the runtime side, `AgentDef` mirrors `ProcessDef` (Bindable/Chainable, has a name, owner script, params), `AgentFactory` mirrors `ProcessFactory`, and `BaseScript` exposes a `protected void agent(String name, Closure body)` method. The body's call-time behavior throws `UnsupportedOperationException("agent execution not yet implemented")` — this plan only delivers the surface; the next plan delivers execution.

**Tech Stack:** Java 17, Groovy 4, ANTLR4 grammar (`ScriptLexer.g4` / `ScriptParser.g4`), Spock for tests, Gradle multi-module build.

---

## Spec reference

This plan implements §2 (DSL surface, partial — directives + I/O + prompt block parsing), §7 module layout (the core/`nf-lang` portions only), and the AST/runtime-model parts of §3. Sections of the spec deferred to follow-up plans:

- §3 agent loop runtime (langchain4j integration)
- §4 tool bridge
- §5 `agent` config scope wiring at runtime
- §6 LLM error handling
- §7 `nf-agent` plugin

Spec doc: `docs/superpowers/specs/2026-05-04-nextflow-llm-agent-design.md`.

---

## File structure

### Files to create

| Path | Responsibility |
|---|---|
| `modules/nf-lang/src/main/java/nextflow/script/ast/AgentNode.java` | AST node for an agent definition. Extends `MethodNode`. Holds `directives`, `inputs`, `outputs`, `prompt` like `ProcessNodeV2`. |
| `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy` | Runtime model. Mirrors `ProcessDef` shape — name, owner, body. `apply()` throws `UnsupportedOperationException` in this plan. |
| `modules/nextflow/src/main/groovy/nextflow/script/AgentFactory.groovy` | Mirrors `ProcessFactory`. Holds session and owner; constructs `AgentDef` from a body closure. |
| `modules/nextflow/src/test/groovy/nextflow/script/AgentDefTest.groovy` | Spock unit test for `AgentDef`. |
| `modules/nextflow/src/test/groovy/nextflow/script/AgentFactoryTest.groovy` | Spock unit test for `AgentFactory`. |
| `modules/nextflow/src/test/groovy/nextflow/script/parser/v2/AgentScriptLoadingTest.groovy` | Spock end-to-end loader test (Task 15). |
| `modules/nf-lang/src/test/groovy/nextflow/script/parser/AgentParserTest.groovy` | Spock test for the grammar/AST builder, modelled on `ScriptAstBuilderTest`. Kept separate to avoid bloating that file. |

### Files to modify

| Path | Change |
|---|---|
| `modules/nf-lang/src/main/antlr/ScriptLexer.g4` | Add `AGENT : 'agent';` token next to `PROCESS`. |
| `modules/nf-lang/src/main/antlr/ScriptParser.g4` | Add `agentDef` to `scriptDeclaration` alternatives; add `agentDef`, `agentBody`, `agentDirectives`, `agentInputs`, `agentOutputs`, `agentPrompt` rules. |
| `modules/nf-lang/src/main/java/nextflow/script/ast/ScriptNode.java` | Add `private final List<AgentNode> agents`, `addAgent`, `getAgents`. |
| `modules/nf-lang/src/main/java/nextflow/script/ast/ScriptVisitor.java` | Add `void visitAgent(AgentNode node);`. |
| `modules/nf-lang/src/main/java/nextflow/script/ast/ScriptVisitorSupport.java` | Default `visitAgent` impl (no-op). |
| `modules/nf-lang/src/main/java/nextflow/script/parser/ScriptAstBuilder.java` | Handle the new `AgentDefAltContext`; add private `agentDef(...)`, `agentDirectives(...)`, `agentInputs(...)`, `agentOutputs(...)`, `agentPrompt(...)`. Add `"agent"` to `SCRIPT_DEF_NAMES`. |
| `modules/nf-lang/src/main/java/nextflow/script/control/ScriptResolveVisitor.java` | Iterate agents alongside processes/workflows when resolving symbols. |
| `modules/nf-lang/src/main/java/nextflow/script/control/VariableScopeVisitor.java` | Same — agents need variable scoping. |
| `modules/nf-lang/src/main/java/nextflow/script/control/CallSiteCollector.java` | Collect call sites inside agents. |
| `modules/nf-lang/src/main/java/nextflow/script/control/ResolveIncludeVisitor.java` | Add agents to the includable definitions list. |
| `modules/nextflow/src/main/groovy/nextflow/script/BaseScript.groovy` | Add `protected void agent(String name, Closure body)` mirroring `process(...)`. |

The `nf-agent` plugin and the langchain4j integration are **not** touched in this plan.

---

## Pre-flight

- [ ] **Step 0.1: Confirm worktree state**

```bash
git -C /Users/pditommaso/Projects/nextflow status -sb
```

Expected: clean or only `docs/` changes from the spec/plan. If unrelated work is staged, stop and surface it before starting.

- [ ] **Step 0.2: Confirm baseline tests pass**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:test --tests 'nextflow.script.parser.ScriptAstBuilderTest' -q
```

Expected: BUILD SUCCESSFUL. If failing, fix or rebase before proceeding — the new tests will be added to the same module.

---

## Task 1: Failing parser test for a minimal agent definition

This test drives the grammar + AST changes. Created first so every later step has a target to satisfy.

**Files:**
- Create: `modules/nf-lang/src/test/groovy/nextflow/script/parser/AgentParserTest.groovy`

- [ ] **Step 1.1: Write the failing test**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
package nextflow.script.parser

import nextflow.script.ast.AgentNode
import spock.lang.Specification
import test.TestUtils

class AgentParserTest extends Specification {

    def setupSpec() {
        TestUtils.beforeSpec()
    }

    def 'should parse a minimal agent definition'() {
        when:
        def script = TestUtils.parse('''\
            nextflow.enable.types = true

            agent eval_agent {
                model 'openai/gpt-5-mini'
                instruction 'You are helpful.'
                tools()
                maxIterations 20

                input:
                    val question

                output:
                    val plan

                prompt:
                """
                Question: ${question}
                """
            }
            '''.stripIndent())

        then:
        script.agents.size() == 1
        def node = script.agents[0] as AgentNode
        node.name == 'eval_agent'
        node.inputs.length == 1
        node.inputs[0].name == 'question'
    }
}
```

- [ ] **Step 1.2: Verify the test fails to compile**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:compileTestGroovy -q
```

Expected: compilation fails — `AgentNode`, `script.agents`, etc. don't exist yet. This confirms the test is exercising new surface.

- [ ] **Step 1.3: Commit the failing test**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nf-lang/src/test/groovy/nextflow/script/parser/AgentParserTest.groovy
git -C /Users/pditommaso/Projects/nextflow commit -s -m "test: add failing parser test for agent definition"
```

---

## Task 2: Add the `AGENT` lexer token

**Files:**
- Modify: `modules/nf-lang/src/main/antlr/ScriptLexer.g4`

- [ ] **Step 2.1: Add the AGENT token**

Find the line containing `PROCESS         : 'process';` (around line 345). Add immediately above it:

```antlr
AGENT           : 'agent';
```

Result (context):

```antlr
// -- agent definition
AGENT           : 'agent';

// -- process definition
PROCESS         : 'process';
```

- [ ] **Step 2.2: Verify the lexer compiles**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:generateGrammarSource -q
```

Expected: BUILD SUCCESSFUL. ANTLR regenerates lexer/parser sources cleanly.

- [ ] **Step 2.3: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nf-lang/src/main/antlr/ScriptLexer.g4
git -C /Users/pditommaso/Projects/nextflow commit -s -m "lang: add AGENT token to script lexer"
```

---

## Task 3: Add `agentDef` parser rules

**Files:**
- Modify: `modules/nf-lang/src/main/antlr/ScriptParser.g4`

- [ ] **Step 3.1: Add the alternative to `scriptDeclaration`**

Locate the alternatives block (around line 117 — `processDef #processDefAlt` etc.) and insert above `processDef`:

```antlr
    |   agentDef                    #agentDefAlt
```

Resulting block:

```antlr
    :   functionDef                  #functionDefAlt
    |   agentDef                    #agentDefAlt
    |   processDef                  #processDefAlt
    |   workflowDef                 #workflowDefAlt
    ...
```

- [ ] **Step 3.2: Add the agent grammar rules**

Append the following rules immediately after the `processStub` rule block (search for `processStub` and add after the last process rule, before the workflow rules around line 357):

```antlr
// -- agent definition
agentDef
    :   AGENT name=identifier nls LBRACE
        body=agentBody?
        sep? RBRACE
    ;

agentBody
    :   (sep agentDirectives)?
        (sep agentInputs)?
        (sep agentOutputs)?
        sep agentPrompt
    ;

agentDirectives
    :   statement (sep statement)*
    ;

agentInputs
    :   INPUT COLON nls processInput (sep processInput)*
    ;

agentOutputs
    :   OUTPUT COLON nls processOutput (sep processOutput)*
    ;

agentPrompt
    :   PROMPT COLON nls statement
    ;
```

The `agentInputs`/`agentOutputs` rules deliberately reuse `processInput` and `processOutput` — agents accept the same typed parameter shape as processV2.

- [ ] **Step 3.3: Add the `PROMPT` lexer token**

In `modules/nf-lang/src/main/antlr/ScriptLexer.g4`, find `INPUT : 'input';` and add nearby (alphabetical placement is fine):

```antlr
PROMPT          : 'prompt';
```

- [ ] **Step 3.4: Verify the parser regenerates**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:generateGrammarSource -q
```

Expected: BUILD SUCCESSFUL. If ANTLR reports rule conflicts, surface and stop — do not paper over.

- [ ] **Step 3.5: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nf-lang/src/main/antlr/ScriptLexer.g4 modules/nf-lang/src/main/antlr/ScriptParser.g4
git -C /Users/pditommaso/Projects/nextflow commit -s -m "lang: add agentDef grammar rules"
```

---

## Task 4: Create the `AgentNode` AST class

**Files:**
- Create: `modules/nf-lang/src/main/java/nextflow/script/ast/AgentNode.java`

- [ ] **Step 4.1: Write `AgentNode.java`**

```java
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
package nextflow.script.ast;

import org.codehaus.groovy.ast.ClassHelper;
import org.codehaus.groovy.ast.MethodNode;
import org.codehaus.groovy.ast.Parameter;
import org.codehaus.groovy.ast.stmt.EmptyStatement;
import org.codehaus.groovy.ast.stmt.Statement;

/**
 * AST node for an agent definition.
 *
 * Mirrors the shape of {@link ProcessNodeV2} but with a simpler body:
 * directives, typed inputs, typed outputs, and a prompt block.
 * No script/exec/stub/when/stage/topic — execution is delegated to the
 * agent runner (see nf-agent plugin).
 */
public class AgentNode extends MethodNode {

    public final Statement directives;
    public final Parameter[] inputs;
    public final Statement outputs;
    public final Statement prompt;

    public AgentNode(String name, Statement directives, Parameter[] inputs, Statement outputs, Statement prompt) {
        super(name, 0, ClassHelper.OBJECT_TYPE, inputs, ClassHelper.EMPTY_TYPE_ARRAY, EmptyStatement.INSTANCE);
        this.directives = directives;
        this.inputs = inputs;
        this.outputs = outputs;
        this.prompt = prompt;
    }
}
```

- [ ] **Step 4.2: Verify it compiles**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:compileJava -q
```

Expected: BUILD SUCCESSFUL.

- [ ] **Step 4.3: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nf-lang/src/main/java/nextflow/script/ast/AgentNode.java
git -C /Users/pditommaso/Projects/nextflow commit -s -m "lang: add AgentNode AST class"
```

---

## Task 5: Extend `ScriptNode` to hold agents

**Files:**
- Modify: `modules/nf-lang/src/main/java/nextflow/script/ast/ScriptNode.java`

- [ ] **Step 5.1: Add the `agents` collection and accessors**

Find the existing `processes` field (search for `private final List<ProcessNode> processes`). Add an analogous field next to it:

```java
private final List<AgentNode> agents = new ArrayList<>();
```

Find the `addProcess` method (around line 150). Add immediately below:

```java
public void addAgent(AgentNode agentNode) {
    agents.add(agentNode);
}
```

Find the `getProcesses` getter (search for `public List<ProcessNode> getProcesses`). Add an analogous getter:

```java
public List<AgentNode> getAgents() {
    return agents;
}
```

- [ ] **Step 5.2: Verify it compiles**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:compileJava -q
```

Expected: BUILD SUCCESSFUL.

- [ ] **Step 5.3: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nf-lang/src/main/java/nextflow/script/ast/ScriptNode.java
git -C /Users/pditommaso/Projects/nextflow commit -s -m "lang: track agent definitions on ScriptNode"
```

---

## Task 6: Extend `ScriptVisitor` and `ScriptVisitorSupport`

**Files:**
- Modify: `modules/nf-lang/src/main/java/nextflow/script/ast/ScriptVisitor.java`
- Modify: `modules/nf-lang/src/main/java/nextflow/script/ast/ScriptVisitorSupport.java`

- [ ] **Step 6.1: Add `visitAgent` to the interface**

In `ScriptVisitor.java`, add immediately above `void visitProcess(ProcessNode node);`:

```java
void visitAgent(AgentNode node);
```

- [ ] **Step 6.2: Add the default no-op to `ScriptVisitorSupport`**

Open `ScriptVisitorSupport.java` and find the existing `visitProcess` default. Add immediately above:

```java
@Override
public void visitAgent(AgentNode node) {}
```

- [ ] **Step 6.3: Verify all visitor implementations still compile**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:compileJava -q
```

Expected: BUILD SUCCESSFUL. If any subclass fails because it doesn't override `visitAgent`, decide per case: most should inherit the default by extending `ScriptVisitorSupport`. If a class implements the interface directly, add a no-op override.

- [ ] **Step 6.4: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nf-lang/src/main/java/nextflow/script/ast/ScriptVisitor.java modules/nf-lang/src/main/java/nextflow/script/ast/ScriptVisitorSupport.java
git -C /Users/pditommaso/Projects/nextflow commit -s -m "lang: add visitAgent to script visitor"
```

---

## Task 7: Implement `agentDef` in `ScriptAstBuilder`

**Files:**
- Modify: `modules/nf-lang/src/main/java/nextflow/script/parser/ScriptAstBuilder.java`

- [ ] **Step 7.1: Add `"agent"` to `SCRIPT_DEF_NAMES`**

Find the constant (around line 211):

```java
private static final List<String> SCRIPT_DEF_NAMES = List.of("process", "workflow", "output");
```

Replace with:

```java
private static final List<String> SCRIPT_DEF_NAMES = List.of("agent", "process", "workflow", "output");
```

- [ ] **Step 7.2: Add the `AgentDefAltContext` branch**

Find the `else if( ctx instanceof ProcessDefAltContext pdac )` branch (around line 319). Add a parallel branch immediately above it:

```java
else if( ctx instanceof AgentDefAltContext adac ) {
    var node = agentDef(adac.agentDef());
    moduleNode.addAgent(node);
}
```

- [ ] **Step 7.3: Add the `agentDef` builder method**

Add the following private methods to `ScriptAstBuilder.java`. Place them immediately after `processStub(...)` for proximity to the related code:

```java
private AgentNode agentDef(AgentDefContext ctx) {
    var name = identifier(ctx.name);
    if( ctx.body == null )
        return invalidAgent("Missing agent body", ctx);
    if( ctx.body.agentPrompt() == null )
        return invalidAgent("Missing `prompt:` section", ctx);

    var directives = agentDirectives(ctx.body.agentDirectives());
    var inputs = agentInputs(ctx.body.agentInputs());
    var outputs = agentOutputs(ctx.body.agentOutputs());
    var prompt = agentPrompt(ctx.body.agentPrompt());

    var node = ast(new AgentNode(name, directives, inputs, outputs, prompt), ctx);
    return node;
}

private AgentNode invalidAgent(String message, AgentDefContext ctx) {
    var empty = EmptyStatement.INSTANCE;
    var result = ast(new AgentNode("", empty, Parameter.EMPTY_ARRAY, empty, empty), ctx);
    collectSyntaxError(new SyntaxException(message, result));
    return result;
}

private Statement agentDirectives(AgentDirectivesContext ctx) {
    if( ctx == null )
        return EmptyStatement.INSTANCE;
    var stmts = ctx.statement().stream()
        .map(this::statement)
        .map(stmt -> checkDirective(stmt, "Invalid agent directive"))
        .collect(Collectors.toList());
    return block(null, stmts);
}

private Parameter[] agentInputs(AgentInputsContext ctx) {
    if( ctx == null )
        return Parameter.EMPTY_ARRAY;
    return ctx.processInput().stream()
        .map(this::processInput)
        .toArray(Parameter[]::new);
}

private Statement agentOutputs(AgentOutputsContext ctx) {
    if( ctx == null )
        return EmptyStatement.INSTANCE;
    var stmts = ctx.processOutput().stream()
        .map(this::processOutput)
        .collect(Collectors.toList());
    return block(null, stmts);
}

private Statement agentPrompt(AgentPromptContext ctx) {
    return statement(ctx.statement());
}
```

All helper methods used here exist in `ScriptAstBuilder.java` with these exact signatures:
- `private Parameter processInput(ProcessInputContext ctx)` (~line 530)
- `private Statement processOutput(ProcessOutputContext ctx)` (~line 648)
- `private Statement statement(StatementContext ctx)` (~line 998)
- `private Statement checkDirective(Statement stmt, String errorMessage)` (~line 721)
- `block(...)`, `ast(...)`, `identifier(...)`, `collectSyntaxError(...)` are widely used in the file.

- [ ] **Step 7.4: Add the import**

Add to the imports section near the existing AST imports:

```java
import nextflow.script.ast.AgentNode;
```

- [ ] **Step 7.5: Verify it compiles and existing tests still pass**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:compileJava -q
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:test --tests 'nextflow.script.parser.ScriptAstBuilderTest' -q
```

Expected: BUILD SUCCESSFUL on both. The existing parser tests cover process/workflow paths — they must remain green.

- [ ] **Step 7.6: Run the new agent parser test**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:test --tests 'nextflow.script.parser.AgentParserTest' -q
```

Expected: PASS.

- [ ] **Step 7.7: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nf-lang/src/main/java/nextflow/script/parser/ScriptAstBuilder.java
git -C /Users/pditommaso/Projects/nextflow commit -s -m "lang: build AgentNode in ScriptAstBuilder"
```

---

## Task 8: Negative parser test — missing `prompt:` section

**Files:**
- Modify: `modules/nf-lang/src/test/groovy/nextflow/script/parser/AgentParserTest.groovy`

- [ ] **Step 8.1: Add the failing test for missing prompt**

Append to `AgentParserTest.groovy`:

```groovy
def 'should report an error for agent without prompt section'() {
    when:
    def errors = TestUtils.check('''\
        nextflow.enable.types = true

        agent broken {
            model 'openai/gpt-5-mini'
            instruction 'x'
            tools()

            input:
                val q
        }
        '''.stripIndent())

    then:
    errors.size() == 1
    errors[0].getOriginalMessage() == 'Missing `prompt:` section'
}
```

(Confirm `TestUtils.check` is the helper that returns a list of `SyntaxException` — see `ScriptAstBuilderTest` for the pattern. If the helper has a different name, adjust.)

- [ ] **Step 8.2: Run it**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:test --tests 'nextflow.script.parser.AgentParserTest' -q
```

Expected: PASS. The error path is already covered by Task 7's `invalidAgent` branch.

- [ ] **Step 8.3: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nf-lang/src/test/groovy/nextflow/script/parser/AgentParserTest.groovy
git -C /Users/pditommaso/Projects/nextflow commit -s -m "test: cover invalid agent definition without prompt"
```

---

## Task 9: Wire agents into resolution and include visitors

These visitors iterate the script-level definitions for symbol resolution and includes. Without these updates, an agent reference in a workflow would resolve as unknown.

**Files:**
- Modify: `modules/nf-lang/src/main/java/nextflow/script/control/ScriptResolveVisitor.java`
- Modify: `modules/nf-lang/src/main/java/nextflow/script/control/VariableScopeVisitor.java`
- Modify: `modules/nf-lang/src/main/java/nextflow/script/control/CallSiteCollector.java`
- Modify: `modules/nf-lang/src/main/java/nextflow/script/control/ResolveIncludeVisitor.java`

- [ ] **Step 9.1: Failing test — agent referenced from a workflow**

In `AgentParserTest.groovy`, add:

```groovy
def 'should resolve an agent reference from a workflow'() {
    when:
    def script = TestUtils.parse('''\
        nextflow.enable.types = true

        agent eval_agent {
            model 'openai/gpt-5-mini'
            instruction 'x'
            tools()

            input:
                val q

            output:
                val r

            prompt:
            """
            ${q}
            """
        }

        workflow {
            channel.of('hi') | eval_agent | view
        }
        '''.stripIndent())

    then:
    script.agents.size() == 1
    script.workflows.size() == 1
    // Workflow body must reference the agent without a "variable not found" error.
    // TestUtils.parse fails if any error is collected.
}
```

- [ ] **Step 9.2: Run it (expected fail)**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:test --tests 'nextflow.script.parser.AgentParserTest' -q
```

Expected: FAIL — `eval_agent` not resolved as a symbol, since visitors only iterate processes/workflows.

- [ ] **Step 9.3: Update `ScriptResolveVisitor`**

Find the loop block (around line 100-103) that already calls `visitWorkflow` then `visitProcess`. Insert an agent loop between them:

```java
for( var workflowNode : sn.getWorkflows() )
    visitWorkflow(workflowNode);
for( var agentNode : sn.getAgents() )
    visitAgent(agentNode);
for( var processNode : sn.getProcesses() )
    visitProcess(processNode);
```

Then add a `visitAgent` method modelled on the existing `visitProcess`. The agent's relevant statements are `node.directives`, `node.outputs`, and `node.prompt` (inputs are `Parameter[]` declared at the method-node level and resolved alongside other parameters). Apply `resolver.transform(...)` / `visit(...)` to each, matching the process pattern.

Add the `import nextflow.script.ast.AgentNode;` at the top of the file.

- [ ] **Step 9.4: Update `VariableScopeVisitor`**

In the loop block around line 113-117, add agents alongside processes (agents extend `MethodNode`, so `declareMethod` accepts them directly):

```java
for( var processNode : sn.getProcesses() )
    declareMethod(processNode);
for( var agentNode : sn.getAgents() )
    declareMethod(agentNode);
```

Add the `import nextflow.script.ast.AgentNode;` import.

- [ ] **Step 9.5: Update `ResolveIncludeVisitor`**

In the `result.addAll(...)` block (around line 174-175), add agents so they can be re-exported through `include { my_agent } from './path.nf'`:

```java
result.addAll(scriptNode.getWorkflows());
result.addAll(scriptNode.getProcesses());
result.addAll(scriptNode.getAgents());
```

`CallSiteCollector` is **not** updated in this plan — it collects compile-time invocations inside workflow bodies, and agents in v1 do not invoke other components from inside their body (no `script:` block). Tool references in `tools(...)` are resolved by `ScriptResolveVisitor` like any other symbol.

- [ ] **Step 9.6: Run the test**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:test --tests 'nextflow.script.parser.AgentParserTest' -q
```

Expected: PASS. Also rerun the full nf-lang test suite to catch regressions:

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:test -q
```

Expected: BUILD SUCCESSFUL.

- [ ] **Step 9.7: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nf-lang/src/main/java/nextflow/script/control/ modules/nf-lang/src/test/groovy/nextflow/script/parser/AgentParserTest.groovy
git -C /Users/pditommaso/Projects/nextflow commit -s -m "lang: include agents in symbol resolution and include visitors"
```

---

## Task 10: Failing runtime test for `AgentDef`

**Files:**
- Create: `modules/nextflow/src/test/groovy/nextflow/script/AgentDefTest.groovy`

- [ ] **Step 10.1: Write the failing test**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
package nextflow.script

import spock.lang.Specification

class AgentDefTest extends Specification {

    def 'should construct an AgentDef with name and body'() {
        given:
        def script = Mock(BaseScript)
        def body = { -> /* placeholder */ } as Closure

        when:
        def agent = new AgentDef(script, body, 'eval_agent')

        then:
        agent.name == 'eval_agent'
        agent.simpleName == 'eval_agent'
        agent.type == 'agent'
    }

    def 'should throw on run() since execution is not yet implemented'() {
        given:
        def script = Mock(BaseScript)
        def body = { -> } as Closure
        def agent = new AgentDef(script, body, 'foo')

        when:
        agent.run(new Object[0])

        then:
        def e = thrown(UnsupportedOperationException)
        e.message.contains('agent execution not yet implemented')
    }

    def 'should clone with a new name'() {
        given:
        def script = Mock(BaseScript)
        def agent = new AgentDef(script, { -> }, 'foo')

        when:
        def renamed = agent.cloneWithName('bar')

        then:
        renamed instanceof AgentDef
        renamed.name == 'bar'
        agent.name == 'foo' // original untouched
    }
}
```

- [ ] **Step 10.2: Verify it fails to compile**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nextflow:compileTestGroovy -q
```

Expected: compilation fails — `AgentDef` does not exist.

- [ ] **Step 10.3: Commit the failing test**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nextflow/src/test/groovy/nextflow/script/AgentDefTest.groovy
git -C /Users/pditommaso/Projects/nextflow commit -s -m "test: add failing AgentDef construction test"
```

---

## Task 11: Implement `AgentDef`

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy`

- [ ] **Step 11.1: Write `AgentDef.groovy`**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
package nextflow.script

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Runtime model for an `agent` definition.
 *
 * Mirrors {@link ProcessDef}'s contract — the script-level binding is callable
 * and chainable in a workflow — but execution is delegated to the future
 * `nf-agent` plugin runner. Calling {@link #apply()} in this plan throws
 * {@link UnsupportedOperationException} until the runner lands.
 */
@Slf4j
@CompileStatic
class AgentDef extends BindableDef implements ChainableDef {

    static final String TYPE = 'agent'

    private BaseScript owner
    private String name
    private String simpleName
    private Closure body

    AgentDef(BaseScript owner, Closure body, String name) {
        this.owner = owner
        this.body = body
        this.name = name
        this.simpleName = name
    }

    @Override
    String getType() { TYPE }

    @Override
    String getName() { name }

    String getSimpleName() { simpleName }

    BaseScript getOwner() { owner }

    Closure getBody() { body }

    @Override
    ComponentDef cloneWithName(String name) {
        def copy = (AgentDef) this.clone()
        copy.@name = name
        return copy
    }

    @Override
    Object run(Object[] args) {
        // Inherited path: BindableDef.invoke_a clones this and calls run().
        throw new UnsupportedOperationException('agent execution not yet implemented')
    }

    Object apply() {
        throw new UnsupportedOperationException('agent execution not yet implemented')
    }
}
```

The class hierarchy mirrors `ProcessDef extends BindableDef implements IterableDef, ChainableDef` — minus `IterableDef` since agents aren't channel iterators in v1. The required abstract method overrides are:
- `getType()` from `ComponentDef` — returns the literal `'agent'` (used in error messages and `toString`).
- `getName()` from `ComponentDef` — returns the agent name.
- `cloneWithName(String)` from `ComponentDef` — clones and renames; standard pattern.
- `run(Object[])` from `BindableDef` — execution entry point. Throws `UnsupportedOperationException` until the runner lands.

`invoke_a` is inherited from `BindableDef` (which calls `run()`); no override needed.

- [ ] **Step 11.2: Run the AgentDef tests**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nextflow:test --tests 'nextflow.script.AgentDefTest' -q
```

Expected: PASS.

- [ ] **Step 11.3: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nextflow/src/main/groovy/nextflow/script/AgentDef.groovy
git -C /Users/pditommaso/Projects/nextflow commit -s -m "core: add AgentDef runtime model (execution stubbed)"
```

---

## Task 12: Failing runtime test for `AgentFactory`

**Files:**
- Create: `modules/nextflow/src/test/groovy/nextflow/script/AgentFactoryTest.groovy`

- [ ] **Step 12.1: Write the failing test**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
package nextflow.script

import nextflow.Session
import spock.lang.Specification

class AgentFactoryTest extends Specification {

    def 'should build an AgentDef from a body closure'() {
        given:
        def session = Mock(Session)
        def script = Mock(BaseScript)
        def factory = new AgentFactory(script, session)

        when:
        def agent = factory.newAgent('eval_agent', { -> })

        then:
        agent instanceof AgentDef
        agent.name == 'eval_agent'
    }
}
```

- [ ] **Step 12.2: Verify it fails to compile**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nextflow:compileTestGroovy -q
```

Expected: compilation fails — `AgentFactory` does not exist.

- [ ] **Step 12.3: Commit the failing test**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nextflow/src/test/groovy/nextflow/script/AgentFactoryTest.groovy
git -C /Users/pditommaso/Projects/nextflow commit -s -m "test: add failing AgentFactory test"
```

---

## Task 13: Implement `AgentFactory`

**Files:**
- Create: `modules/nextflow/src/main/groovy/nextflow/script/AgentFactory.groovy`

- [ ] **Step 13.1: Write `AgentFactory.groovy`**

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
package nextflow.script

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.Session

/**
 * Factory for {@link AgentDef} instances. Counterpart to {@link ProcessFactory}.
 *
 * Deliberately thin in this plan — the heavy lifting (chat-model construction,
 * tool registration, agent loop) lives in the future nf-agent plugin runner.
 */
@Slf4j
@CompileStatic
class AgentFactory {

    private Session session
    private BaseScript owner

    AgentFactory(BaseScript ownerScript, Session session) {
        this.owner = ownerScript
        this.session = session
    }

    AgentDef newAgent(String name, Closure body) {
        new AgentDef(owner, body, name)
    }
}
```

- [ ] **Step 13.2: Run the AgentFactory tests**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nextflow:test --tests 'nextflow.script.AgentFactoryTest' -q
```

Expected: PASS.

- [ ] **Step 13.3: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nextflow/src/main/groovy/nextflow/script/AgentFactory.groovy
git -C /Users/pditommaso/Projects/nextflow commit -s -m "core: add AgentFactory"
```

---

## Task 14: Add `agent(...)` to `BaseScript`

**Files:**
- Modify: `modules/nextflow/src/main/groovy/nextflow/script/BaseScript.groovy`

- [ ] **Step 14.1: Add the `agent` method**

Find the existing `protected void process(...)` method (around line 137). Add the following method immediately above it:

```groovy
/**
 * Define an agent. Mirrors {@link #process(String, Closure)} but constructs an
 * {@link AgentDef} via {@link AgentFactory}.
 */
protected void agent(String name, Closure body) {
    final factory = new AgentFactory(this, session)
    final agent = factory.newAgent(name, body)
    meta.addDefinition(agent)
}
```

`ScriptMeta.addDefinition(ComponentDef)` accepts our `AgentDef` because `BindableDef extends ComponentDef`, so this should compile without further changes.

- [ ] **Step 14.2: Verify the module compiles**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nextflow:compileGroovy -q
```

Expected: BUILD SUCCESSFUL. If `meta.addDefinition` complains about the type, fix `AgentDef`'s supertype as described above and re-run.

- [ ] **Step 14.3: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nextflow/src/main/groovy/nextflow/script/BaseScript.groovy
git -C /Users/pditommaso/Projects/nextflow commit -s -m "core: expose agent { ... } in BaseScript"
```

---

## Task 15: End-to-end smoke test — load a script with an agent

**Files:**
- Create: `modules/nextflow/src/test/groovy/nextflow/script/AgentScriptLoadingTest.groovy`

- [ ] **Step 15.1: Write the smoke test**

Pattern is taken directly from `modules/nextflow/src/test/groovy/nextflow/script/parser/v2/ScriptLoaderV2Test.groovy` ("should run a file script") so the loader entry point matches the v2 parser.

```groovy
/*
 * Copyright 2013-2026, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 */
package nextflow.script.parser.v2

import java.nio.file.Files

import nextflow.Session
import nextflow.script.AgentDef
import nextflow.script.ScriptMeta
import test.Dsl2Spec

class AgentScriptLoadingTest extends Dsl2Spec {

    def 'should load a script with an agent definition'() {
        given:
        def session = new Session()
        def parser = new ScriptLoaderV2(session)
        def file = Files.createTempDirectory('test').resolve('main.nf')
        file.text = '''
            nextflow.enable.types = true

            agent hello_agent {
                model 'openai/gpt-5-mini'
                instruction 'be helpful'
                tools()

                input:
                    val q

                output:
                    val r

                prompt:
                """
                ${q}
                """
            }

            workflow {
                /* no invocation — execution not yet implemented */
            }
            '''.stripIndent()

        when:
        parser.parse(file)

        then:
        def definitions = ScriptMeta.get(parser.script).getDefinitions()
        definitions.any { it instanceof AgentDef && it.name == 'hello_agent' }

        cleanup:
        file.parent.deleteDir()
    }
}
```

- [ ] **Step 15.2: Run it**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nextflow:test --tests 'nextflow.script.AgentScriptLoadingTest' -q
```

Expected: PASS. If it fails, iterate on the failure — typical issues: missing nextflow.enable.types preview gate, missing addDefinition wiring, name collision with builtin.

- [ ] **Step 15.3: Commit**

```bash
git -C /Users/pditommaso/Projects/nextflow add modules/nextflow/src/test/groovy/nextflow/script/AgentScriptLoadingTest.groovy
git -C /Users/pditommaso/Projects/nextflow commit -s -m "test: end-to-end loading of a script with an agent definition"
```

---

## Task 16: Final verification

- [ ] **Step 16.1: Full test suite for both touched modules**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew :nf-lang:test :nextflow:test -q
```

Expected: BUILD SUCCESSFUL. No regressions.

- [ ] **Step 16.2: Quick sanity check — assemble**

```bash
cd /Users/pditommaso/Projects/nextflow && ./gradlew assemble -q
```

Expected: BUILD SUCCESSFUL.

- [ ] **Step 16.3: Confirm git log**

```bash
git -C /Users/pditommaso/Projects/nextflow log --oneline -20
```

Expected: 15 commits from Task 1 → Task 15, all signed off.

---

## What this plan does NOT deliver

These belong to the next plan(s):

1. **Agent execution** — `AgentDef.apply()` and `invoke_a()` throw `UnsupportedOperationException`. A workflow that invokes an agent will fail at runtime with a clear message.
2. **`nf-agent` plugin** — no langchain4j dependency, no chat-model construction, no tool dispatcher.
3. **Tool bridge** — `ModuleSpec` → `ToolSpecification` mapping is not implemented.
4. **`agent` config scope** — no `nextflow.config` `agent { ... }` parsing yet.
5. **Documentation pages** — `docs/reference/process.md` analogue for agents will follow once execution exists; documenting a stub is misleading.

Wiring these up is straightforward once this plan lands: the runner SPI is the single seam where execution plugs in.
