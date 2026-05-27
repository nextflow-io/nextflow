# Type system

- Authors: Ben Sherman
- Status: accepted
- Deciders: Ben Sherman, Paolo Di Tommaso
- Date: 2026-05-1
- Tags: lang, static-types

## Summary

Implement a type system for the Nextflow language: standard types, type annotations, and static type checking.

## Problem Statement

Nextflow is a dynamically-typed language, which means that the types of values are determined at runtime rather than compile-time. This approach allowed Nextflow to iterate rapidly early in its history, and it allowed non-technical users to write pipelines without learning advanced programming techniques.

However, as Nextflow has matured and gained industry adoption, there is an increasing need to provide a first-class experience that supports the development of large production pipelines. Dynamic typing has several limitations in this regard:

- **Poor data modeling.** Users need to be able to *model their domain*, that is, to describe the structure of their data as it flows through a pipeline. However, this cannot be done if the language does not have syntax for expressing types. Users can add comments, but comments cannot be verified by the compiler, so they can quickly become out of sync with the code.

- **Poor error checking.** Users need early and actionable error feedback when developing a Nextflow pipeline. However, a dynamically-typed language tends to report more errors at runtime, which are harder and more expensive to fix than compile-time errors.

The strict syntax parser allows Nextflow to surface many errors at compile-time that it could not do before. However, because Nextflow is dynamically-typed, there are entire categories of errors which cannot be detected at compile-time because the compiler does not have sufficient information.

Therefore, a type system is needed to allow users to model their data domain in a way that can be verified at compile-time. Static typing will make Nextflow code (1) easier for users to read and understand and (2) easier for the compiler to surface errors as early as possible.

## Goals

- Define a standard library of namespaces and types for Nextflow, rather than implicitly delegating to Groovy.

- Introduce type annotations as a first-class feature in the Nextflow language.

- Implement static type checking in the compiler and language server.

## Non-goals

- Removing support for dynamically-typed code -- existing code should continue to work, and it should be possible to use dynamically-typed scripts with statically-typed scripts.

- Replacing the existing runtime types -- the type system should use existing types as much as possible, so that existing code continues to work.

- Implementing null analysis -- while a syntax for nullable types may be introduced, static null analysis will be explored in future efforts.

## Decision

Define the standard types for Nextflow as a focused subset of Groovy types with Nextflow-specific additions. Add support for type annotations (`<name>: <type>`) at every level of the Nextflow language. Use the `?` suffix to denote nullable types (no null analysis yet). Implement static type checking with local variable type inference (lenient by default, optional strict mode). Implement type checking in the language server first, then in the runtime once it is mature.

## Core Capabilities

### Standard library

Nextflow inherits a rich set of types from Groovy (and Java by extension). However, Groovy provides significantly more types and methods than are necessary for writing Nextflow pipelines.

Just as the [strict parser](./20250508-strict-syntax-parser.md) allows us to define Nextflow's syntax independently from Groovy, a standard library is needed to define the set of constants, functions, and types in Nextflow. This way, the user does not need to consult external documentation (e.g. Groovy/Java APIs) to know what they can use in Nextflow.

The standard types are documented [here](https://docs.seqera.io/nextflow/reference/stdlib-types). They are designed based on the following considerations:

- Use existing Java/Groovy types as much as possible: `Boolean`, `Float`, `Integer`, `Path`, `String`, `List`, `Map`, `Set`. These types have been battle-tested by the Java ecosystem and are already used in existing code.

- Provide a minimal and focused set of types. Use reference types over primitive types (e.g. `Integer` over `int`). Use interfaces over implementations (e.g. `List` over `ArrayList`). For each type, include only the methods that are needed. Every additional type and method adds cognitive overhead for users, so they should only be included when they add value.

- Types should be easy to hash (e.g. MD5) so that they can be used with task caching.

- Types should be easy to serialize (e.g. to/from JSON) so that they can be integrated seamlessly with external data storage.

- Encourage composition over inheritance for data modeling. Inheritance is an advanced programming technique that adds a lot of complexity to the type checker and isn't needed in Nextflow anyway. When it is useful, the standard library uses *traits* to model shared behaviors among multiple concrete types. For example, the `Iterable` trait is used by the collection types (`List`, `Set`, `Bag`).

Nextflow types can be parameterized like Java/Groovy, e.g. `Channel<Record>` or `Map<String,Integer>`. Wildcards are supported (e.g. `List<?>`), but lower/upper bounds are not (e.g. `List<? extends Record>`). Since inheritance is not used, lower/upper bounds are not needed, which simplifies both the language and the type checker.

The standard library also defines a set of *namespaces*, which provide standalone constants and functions. They are documented [here](https://docs.seqera.io/nextflow/reference/stdlib-namespaces). These namespaces can be extended in the future as needed.

Namespaces should be used instead of static methods. For example, channel factories should be called using the `channel` namespace instead of the `Channel` type.

### User-defined types

Users can define their own types using *enum types* and *record types*.

An enum type models a choice between a set of options:

```groovy
enum Color {
    RED,
    GREEN,
    BLUE
}
```

A record type models a composition of multiple values:

```groovy
record FastqPair {
    id: String
    fastq_1: Path
    fastq_2: Path
}
```

The design and implementation of record types is described in the [Record types ADR](./20260306-record-types.md).

User-defined types can be included across modules like other script definitions, and they can be used in type annotations like the standard types:

```groovy
include { Color ; FastqPair } from './types.nf'

workflow hello {
    take:
    samples: Channel<FastqPair>
    color: Color

    // ...
}
```

Users can model arbitrary data structures by composing enum types, record types, and the standard types.

### Type annotations

Type annotations are needed to express the types of values at every level of the Nextflow language.

For example:

```groovy
// pipeline params and outputs
params {
    ids: List<String>
    index: Path
}

workflow {
    // ...
}

output {
    samples: Channel<Sample> {
        // ...
    }

    multiqc_report: Path {
        // ...
    }
}

// workflow takes and emits
workflow SRA {
    take:
    ids: Channel<String>

    main:
    // ...

    emit:
    samples: Channel<Sample>
}

// process inputs and outputs
process FASTQC {
    input:
    reads: List<Path>

    output:
    logs: Path = file('fastqc_logs')

    // ...
}

// function parameters and returns
def isSraId(id: String) -> Boolean {
    // local variable
    def x: String = 'hello'
}
```

The `<name>: <type>` syntax emphasizes the name as the identity and the type as an optional enhancement. It is increasingly used by modern languages such as Python, Rust, and TypeScript, so it will be familiar to most users coming from other languages.

Notably, Java and Groovy use `<type> <name>`. This syntax was never formally part of the Nextflow language, but some advanced users do use them for functions and local variables. To ease the transition for these users, the strict parser should be able to recognize Groovy type annotations and convert them to Nextflow type annotations.

Any type annotation can be marked as *nullable* by appending a `?`:

```groovy
def x_opt: String? = null
```

Static null checking will not be implemented yet, but it can be used for documentation purposes. Null checking may be implemented at runtime, such as the `params` block and typed process inputs.

The implementations of type annotations for specific language features are described in the following ADRs:

- [Workflow params](./20250825-workflow-params.md)
- [Workflow outputs](./20251020-workflow-outputs.md)
- [Typed processes](./20251017-typed-processes.md)
- [Typed workflows](./20260310-typed-workflows.md)

### Static type checking

Static type checking is an additional compilation phase in which the type of every expression is resolved and validated. It is enabled by having a standard set of types and syntax for using those types in the language. Type checking can catch many errors that would otherwise be caught at runtime, such as a function being called with the wrong number of arguments.

Type checking should provide an additional layer of validation without imposing undue burden on the user.

- The type checker should infer the types of expressions where possible, in order to minimize the amount of extra work required to benefit from type checking. Users should only need to specify the types of inputs -- pipeline parameters, workflow takes, process inputs, etc -- since the type checker can infer all downstream code from there.

- When the type of an expression cannot be resolved, the type checker should not immediately treat it as an error. Depending on the context, it may be more appropriate to ignore it or report a warning instead. This way, users can adopt static typing progressively rather than taking an "all-or-nothing" approach.

Type checking will be implemented using a similar approach as the strict parser:

1. Introduce type checking in the language server as a non-blocking extra layer of validation
2. Introduce type checking in the linter (`nextflow lint`) to enable CLI usage
3. Introduce type checking in the runtime (`nextflow run`)

This phased approach allows us to refine the type checker based on real usage before we enable it across the board.

We may want to give users some control over the strictness of the type checker. For example, an optional "strict" type checking mode could provide maximum validation (treat all ambiguous types as errors), whereas the default mode would take a balanced approach (flag ambiguous types and other non-blocking issues as warnings).

Examples of type checking:

```groovy
// binary expressions
2 + '2'                     // error: `+` operator is not defined for Integer and String
2 + 2                       // ok: Integer + Integer → Integer

// property accesses
workflow.output_dir         // error: `workflow` namespace has no property called `output_dir`
workflow.outputDir          // ok: workflow.outputDir -> Path

// method calls
'hello'.toUpper()           // error: String has no method called `toUpper`
'hello'.toUpperCase()       // ok: String::toUpperCase() -> String

// call arguments
env()                       // error: `env` function expects 1 argument but received 0
env(42)                     // error: Integer argument is not compatible with String parameter
env('HELLO')                // ok: env(String) -> String

// variable declarations
def x: Integer = 'hello'    // error: Integer `x` cannot be assigned to a String
def y: Integer = 42         // ok

// higher order functions
channel.of(1, 2, 3, 4, 5)   // Channel<Integer>
    .map { v -> v * v }     // Channel<T> map((T) -> R) -> Channel<R>
    .view()                 // Channel<Integer>

channel.of('a', 'b', 'c')   // Channel<String>
    .map { v -> v * v }     // error: `*` operator is not defined for String and String
```

### Static typing and the runtime

The type system exists as a separate layer from the runtime. That is, the type system defines the standard namespaces and types as interfaces, and the runtime types need only be compatible with those interfaces. Nextflow-specific types (e.g. `Duration`, `MemoryUnit`) can implement the corresponding type interface directly. Types inherited from Java (e.g. `Integer`, `String`, `List`) are modeled in the type system using "shim" interfaces, which are used during type checking and discarded prior at runtime.

This approach allows us to introduce type checking without modifying runtime behavior. Nextflow still generates dynamically-typed Groovy code, but this code is validated by the type checker, so it is effectively statically-typed. This approach is similar to [TypeScript](https://www.typescriptlang.org/), where TypeScript code is type-checked at compile-time and converted to JavaScript at runtime by discarding type annotations. Static compilation (e.g. `@CompileStatic`) may be explored in the future, but it is not clear whether the performance improvement of statically-compiled code is worth the extra compilation time and additional constraints on the generated Groovy code.

### Migration plan

Static typing is a significant change to the way that Nextflow pipelines are written. While the *type system* itself is largely a formalization of existing types, *static type checking* requires pipeline code to be more verbose and precise. As a result, it is important that we ease the migration for existing users as much as possible.

To that end, static typing is an opt-in feature that is enabled per-script using the `nextflow.enable.types` feature flag. This flag allows the use of typed processes and typed workflows, and it enables static type checking at compile-time. Typed modules and legacy modules can be used in the same pipeline.

This approach has several benefits:

- Dynamically-typed code can still be used, indefinitely.

- Existing code can be migrated to static typing with minimal changes. New features such as records and workflow outputs can be adopted independently of static typing.

- Existing code can be migrated one module at a time, since typed modules and legacy modules can be used together.

While we hope that the vast majority of users adopt static typing in the near future, we expect that legacy code will remain for a long time. The above approach will allow us to support both coding styles for the long-term as first-class citizens of the Nextflow language.

The feature flag and interoperability rules are described in the [Typed workflows ADR](./20260310-typed-workflows.md).

### Java and Groovy classes

Nextflow allows the use of arbitrary Java/Groovy classes as an escape hatch for gaps in the standard library. However, this pattern is problematic because it is based on the Nextflow runtime classpath, which is difficult to validate at compile-time and can change across releases.

This pattern should eventually be deprecated and disallowed. Common use cases such as reading/writing JSON and HTTP requests should be supported by the standard library instead.

Beyond that, users can define arbitrary Groovy code in the `lib` directory or in a plugin. With a plugin, the user can specify their dependencies explicitly rather than relying on the Nextflow runtime classpath.

## Alternatives

### Optional type

There are multiple ways to deal with nullity in a programming language:

- Disallow `null` and provide an explicit "optional" type (e.g. `Option` in Rust, `Maybe` in Haskell)

- Allow `null` and provide a way to mark values as nullable (e.g. `String?` in Swift, `prop?: string` in TypeScript)

Java 8 introduced an `Optional` type, but it has limited utility because Java already has `null` and the standard library is already designed around this fact. Instead, many Java libraries mark nullable fields and methods with the `@Nullable` annotation, which language servers can use to perform null checking. This annotation is an informal equivalent to the `?` suffix in TypeScript. In fact, Java will likely introduce an explicit `?` suffix for this purpose as part of [Project Valhalla](https://openjdk.org/projects/valhalla/).

While an explicit optional type is easier to formally verify, a `?` suffix is much easier to use. Consider the difference between the following hypothetical examples in Nextflow.

Declaring a nullable variable:

```groovy
// Optional
def name: Optional<String> = Optional.empty()

// ? suffix
def name: String? = null
```

Safe navigation:

```groovy
// Optional
def city: Optional<String> = Optional.ofNullable(user)
    .map { u -> u.address }
    .map { a -> a.city }
def display: String = city.orElse('unknown')

// ? suffix
def city: String? = user?.address?.city
def display: String = city ?: 'unknown'
```

Conditional execution:

```groovy
// Optional
def email_opt: Optional<String> = Optional.ofNullable(user.email)
email_opt.ifPresent { email ->
    sendEmail(email)
}

// ? suffix
def email: String? = user.email
if( email ) {
    sendEmail(email)
}
```

In all cases, the `Optional` type makes the code more verbose and unnatural, whereas the `?` suffix doesn't change the code at all. The `?` suffix also builds naturally on the `?.` and `?:` operators from Groovy.

The `?` suffix is more difficult to validate at compile-time, but the superior developer experience is well worth the cost. This trade-off is the reason why null analysis is not in scope for this ADR.

## Links

- [Nextflow gotchas](https://midnighter.github.io/nextflow-gotchas/)
