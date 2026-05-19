# Strict syntax parser

- Authors: Ben Sherman
- Status: accepted
- Deciders: Ben Sherman, Paolo Di Tommaso
- Date: 2025-05-08
- Tags: lang, parser

## Summary

Implement a formal grammar and parser for Nextflow scripts and config files.

## Problem Statement

Nextflow was originally developed as a Groovy DSL. Nextflow scripts and config files are parsed and executed as Groovy scripts with AST transforms for Nextflow-specific concepts such as process definitions, include statements, and config blocks.

This approach has several structural limitations:

- **No formal grammar.** Because scripts are parsed as arbitrary Groovy, there is no canonical definition of the Nextflow language. Any Groovy syntax is accepted, even those that are unsafe or nonsensical in a pipeline context.

- **Poor error messages.** Syntax and semantic errors are surfaced as Groovy compiler errors, which do not understand Nextflow concepts such as processes and workflows. Users receive messages that reference internal Groovy machinery rather than familiar Nextflow concepts.

- **Poor developer tooling.** A language server, formatter, or linter cannot operate accurately on a language defined only by a series of post-hoc AST rewrites. For example, there is no clear AST node type for a process or workflow -- they are represented as opaque method calls.

- **Language evolution is coupled to Groovy.** Every new Nextflow feature must be expressed as valid Groovy syntax, even though the best way to express a concept in Nextflow might not be valid in Groovy. This hinders our ability to evolve the Nextflow language based on the needs of Nextflow users.

## Goals

- Define a formal grammar for the Nextflow language that is independent of Groovy.

- Implement a Nextflow parser that can be used both by the runtime and developer tooling (linter, formatter, language server).

- Provide actionable, Nextflow-aware error messages for syntax and semantic errors.

- Preserve backwards compatibility with existing code as much as possible (allow breaking changes for unsupported Groovy syntax).

## Non-goals

- Replacing the Groovy runtime -- Nextflow code should still compile and run as Groovy scripts.

- Removing the legacy parser immediately -- both parsers are supported during the transition period.

- Introducing new language features -- the first iteration of the strict parser should focus on supporting existing functionality, with new features addressed in future efforts.

## Decision

Implement a parser for Nextflow code, using an ANTLR grammar based on the Groovy 4 grammar. The parser should check for syntax and semantic errors and produce the same Groovy AST as the legacy parser. The parser should be implemented in a shared module (`nf-lang`) that can be used by the runtime, linter, formatter, and language server. The runtime should be able to use either the strict parser or legacy parser based on an environment variable (`NXF_SYNTAX_PARSER=v1|v2`).

## Core Capabilities

### Module boundaries

The parser is implemented in the `nf-lang` module, which is independent of the Nextflow runtime. This way, it can be used by the linter, language server, and runtime -- the language server can use `nf-lang` without importing the entire runtime.

### Parsing

The Nextflow grammar (`ScriptParser.g4`) is derived from the Groovy grammar, with additional productions for Nextflow-specific concepts such as process and workflow definitions. The grammar is defined in ANTLR and compiled to Java at build time.

The AST builder (`ScriptAstBuilder`) parses and constructs the abstract syntax tree (AST) for a Nextflow script. The AST at this stage is a Groovy AST (`ASTNode`) with several Nextflow-specific node types:

- `ScriptNode`: top-level scripts
- `ProcessNode`: process definition
- `WorkflowNode`: workflow definition
- `IncludeNode`: `include` declaration

These nodes extend standard Groovy AST types so that they can be used seamlessly with existing Groovy infrastructure during AST analysis.

### Semantic analysis

AST construction is followed by *semantic analysis*, in which the AST is inspected for various kinds of errors. Semantic analysis consists of multiple phases:

1. **Include resolution** (`ResolveIncludeVisitor`) resolves each `include` declaration, validating that the included module exists and the included names exist in the module.

2. **Name checking** (`ScriptResolveVisitor`) resolves variable references, function calls, and type names, reporting errors for any undefined references.

3. **Type checking** (`TypeCheckingVisitor`) infers and validates the types of expressions. *NOTE: The initial implementation only performs minimal checking on process and workflow calls. Full type checking will be addressed in future efforts.*

Errors are collected and reported with specific locations and messages.

### Module resolution

The legacy parser resolves include declarations at runtime, compiling each new module when the including script is executed.

The strict parser, however, resolves all included modules eagerly during compilation, so that it can resolve references across modules. When the main script is compiled, it recursively resolves and parses all included modules (`ModuleResolver`) prior to semantic analysis.

### Code generation

Once a script passes semantic analysis without any errors, it is converted to a "pure" Groovy AST (`ScriptToGroovyVisitor`). For example, process nodes are converted into method calls (`BaseScript::process()`) that register the process in the underlying runtime session (`Session`).

The resulting Groovy AST is identical to the one produced by the legacy parser. The validity of this AST is guaranteed by the layers of error checking in the strict parser. All subsequent Groovy compilation and execution steps proceed unchanged.

Each script is compiled to a Groovy class with a unique name (e.g. `_nf_script_abc123`) derived from its URI, in order to avoid class name collisions.

### Config parsing

A separate parser is implemented for config files, including a grammar (`ConfigParser.g4`) and AST builder (`ConfigAstBuilder`). It supports the same basic syntax for expressions and statements, but uses different top-level primitives (config assignments, config blocks, config includes).

Semantic analysis for config files also involves validating config options, for example:

```groovy
process.cpu = 8   // error: invalid config option `process.cpu`
process.cpus = 8  // ok
```

This validation requires config options to be defined in a way that can be extracted at (Nextflow) build-time and provided to the Nextflow compiler and language server.

Whereas the legacy parser evaluates config directly, the strict parser compiles config statements into an internal DSL. For example, given the following snippet:

```groovy
foo.bar = 'hello'

foo {
  bar = 'hello'
}

includeConfig 'foobar.config'
```

The strict parser produces the following Groovy code:

```groovy
assign(['foo', 'bar'], 'hello')

block('foo', {
  assign(['bar'], 'hello')
})

includeConfig('foobar.config')
```

This transformation gives the runtime (`ConfigDsl`) much more control when evaluating config statements, making it easier to maintain and less prone to bugs.

### Developer tooling

The linter (`nextflow lint`) uses the same components described above to perform parsing and semantic analysis. It does not emit anything beyond Nextflow AST because it does not need to execute the code.

The language server takes a similar approach, with a focus on caching and live editing. Third-party tools can similarly use `nf-lang` to parse, analyze, and manipulate Nextflow code.

### Backwards compatibility

The Nextflow runtime can use the v1 parser (legacy) or v2 parser (strict) based on an environment variable:

```bash
# default in Nextflow 25.04
export NXF_SYNTAX_PARSER=v1

# default in Nextflow 26.04
export NXF_SYNTAX_PARSER=v2
```

In principle, the strict parser is just a strict implementation of the Nextflow language, so existing Nextflow code should continue to work with it.

In practice, many users have adopted patterns that are allowed in Groovy but are not intended for use in Nextflow (e.g. classes, annotations). Some syntax patterns are technically valid but redundant (e.g. lambdas vs closures). Overall, there is a great deal of ambiguity in what is or isn't considered part of the Nextflow language.

To that end, we have defined a language specification for Nextflow which enumerates every syntax construct that is supported, as well as a guide that describes the most common patterns that must be rewritten to comply with the strict parser.

In most cases, these changes are minor and easy to fix. As a last resort, any non-compliant code can be moved to the `lib` directory or a plugin.

### Future language features

A custom parser allows us to design new language features using any syntax, as long as it can be compiled to Groovy AST. New language features should be implemented only in the strict parser to reduce maintenance burden and encourage users to migrate.

## Alternatives

### Improve legacy parser

We could commit to the legacy parser and focus on improving it:

- Check for common mistakes and invalid Groovy syntax during AST analysis

- Build tooling (linter, language server) around the legacy parser and rely on it to improve the developer experience

This approach is fundamentally limited because we have no control over parsing or AST construction, only AST analysis. Some errors are much harder, or even impossible, to capture during AST analysis vs parsing. Furthermore, we would always be restricted to Groovy syntax when considering new language features, which makes it difficult to design new features for our needs.

A custom parser is more work up front, but will be easier to maintain over time and gives us significantly more room to extend the language in the future. We also minimize the potential maintenance burden by continuing to use the rest of the Groovy runtime -- Nextflow only handles the transformation from source code to Groovy AST.

### Replace DSL with an SDK

Maintaining a custom language is a lot of work. DSLs are a low-effort alternative, but often have a poor developer experience as a consequence (poor error messages, poor developer tooling, etc). We could abandon the idea of Nextflow as a language and provide it as an SDK instead.

For example, consider the following Nextflow script:

```groovy
process sayHello {
    cpus 2
    memory 4.GB

    input:
    val greeting

    output:
    stdout

    script:
    """
    echo '${greeting} world!'
    """
}

workflow {
    ch_greetings = channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
    sayHello(ch_greetings).view()
}
```

With a Nextflow SDK, this example might be written as follows:

```groovy
import nextflow.Nf

@Nf.Process(
    cpus = 2,
    memory = '4 GB'
)
def echo(String greeting) {
    return [
        """
        echo '${greeting} world!'
        """,
        Nf.stdout()
    ]
}

@Nf.Workflow
def main() {
    ch_greetings = channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
    sayHello(ch_greetings).view()
}
```

A Nextflow SDK has many potential advantages:

- SDKs have a much more "normal" developer experience than DSLs -- errors can still happen, but they are the same "kinds" of errors encountered with any other library.

- SDKs are just libraries in an existing language, so they can be used with existing developer tools for that language.

- SDks can be exposed to multiple programming languages via language bindings, so Nextflow could be used directly in the user's preferred language (e.g. Python).

The problem is that a Nextflow SDK would be significantly more verbose. The above example is extremely simple, but quickly becomes large and unreadable as the pipeline grows. This trade-off may be acceptable to experienced software engineers, but Nextflow is specifically designed for "domain experts" -- non-technical users who benefit from having a "domain-specific" language.

Users who prefer SDKs already have several good options for workflow systems, such as Airflow, Prefect, and Dagster. Rewriting Nextflow as an SDK is unlikely to attract these users. In fact, it would likely only alienate Nextflow's existing user base.

A custom workflow language is harder to develop and maintain, but it is essential to Nextflow's identity. Pipelines written with Nextflow are much more concise, and many of Nextflow's downsides can be addressed with better documentation and developer tooling.

### Replace Groovy DSL with a Python DSL

Nextflow is a custom workflow language based on Groovy. But most users use Python in their regular work, and are not familiar at all with Groovy. Nextflow might attract many more users if it were based on Python instead of Groovy.

This question has multiple aspects that are worth unpacking:

- **Syntax and semantics.** Groovy and Python use different syntax and names to express the same concepts (curly braces vs indentation, `'hello'.length()` vs `len('hello')`). If Nextflow syntax was based on Python, users would be able to learn Nextflow more easily.

- **Runtime interoperability.** Nextflow code runs on the Groovy runtime, so it can't call Python code natively. It can only call Python scripts through a process, which requires extra ceremony to pass data between Nextflow and Python. If Nextflow code ran on Python, users could use Python libraries directly in Nextflow code.

It may be possible to define a Python-based DSL for Nextflow. For example, the aforementioned "Hello World" example might be written as follows in a Python-based DSL:

```python
process sayHello:
  directives:
    cpus = 2
    memory = 4.GB

  input:
    greeting: str

  output:
    stdout()

  script:
    f'''
    echo '{greeting} world!'
    '''

workflow:
  ch_greetings = channel.of('Bonjour', 'Ciao', 'Hello', 'Hola')
  sayHello(ch_greetings).view()
```

This code could be compiled to pure Python and executed with [GraalPy](https://www.graalvm.org/python/). This way, Nextflow could provide a language that uses Python syntax and allows Python imports, but everything would continue to run on the JVM.

However, such a change comes with significant challenges:

- Some Nextflow idioms may not fit well in a Python-based syntax. For example, channel operators make heavy use of closures, which can span multiple lines, but Python lambdas can only specify a single expression. There are workarounds, but the result may be less readable than before.

- A great deal of Nextflow code has already been written and tested in the current syntax. Introducing a new syntax would require either the entire community to migrate, or, more likely, both syntaxes to be supported long-term. Either way introduces a lot of pain and/or confusion for the community.

- Nextflow would still be a *Python-based DSL*, not pure Python. It would still need its own parser and its own developer tooling, just one that targets Python AST instead of Groovy AST.

- Nextflow uses a functional-reactive programming model (a.k.a. channels and operators) for which there is no native equivalent in Python. This is another reason why users find Nextflow unfamiliar, and it won't go away by simply replacing Groovy syntax with Python syntax.

In summary, a Python DSL for Nextflow is an extreme solution that only addresses a small aspect of the developer experience. Since Nextflow is a custom language, the developer experience ultimately depends on the quality of the parser and developer tooling, regardless of whether the syntax is based on Groovy or Python.

### Use a static configuration language

The Nextflow configuration language has a simple and declarative syntax that supports dynamic expressions where appropriate. However, the legacy parser evaluates config files as Groovy scripts, allowing arbitrary scripting constructs such as if statements and functions. As a result, Nextflow config can become very complicated in practice, leading to a poor developer experience.

An alternative would be to use a static language such as JSON, YAML, or TOML. These languages are ideal for configuration because they are simple and predictable. Nextflow idioms such as includes and profiles could be implemented with special conventions. This would make config files easier to debug.

However, an essential aspect of Nextflow config is the ability to use dynamic expressions, such as for a process directive:

```groovy
process.memory = { 8.GB * task.attempt }
```

This feature would be difficult to replicate in a static language like YAML. Dynamic expressions could be specified as plain-text, but they would not benefit from Nextflow-specific tooling (syntax highlighting, linting).

Much of the difficulty that users experience with Nextflow config comes from (1) lack of developer tooling and (2) gaps in native functionality. The Nextflow config language is simple, but a strict parser is needed to enforce this simplicity and provide clear error feedback. Features such as the `resourceLimits` directive and workflow outputs can eliminate prior ad-hoc solutions that involved a lot of complex Nextflow config.

## Links

- Community discussion: [#3107](https://github.com/nextflow-io/nextflow/discussions/3107)
- Implementation: [#4613](https://github.com/nextflow-io/nextflow/pull/4613), [#4744](https://github.com/nextflow-io/nextflow/pull/4744)
- [Nextflow language specification](../docs/reference/syntax.md)
- [Preparing for strict syntax](../docs/strict-syntax.md)
- Refined by: [Plugin Spec](./20250922-plugin-spec.md)
