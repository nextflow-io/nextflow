# Preserve comments in formatter

- Authors: Ben Sherman, Phil Ewels
- Status: accepted
- Date: 2026-08-05
- Tags: lang, formatter, parser

## Summary

The formatter can silently delete, move, or corrupt comments, because comments are only captured for a hand-picked set of AST node positions. Replace this allowlist with an algorithm that attaches *every* comment in the source to an AST node, ensuring that every comment is preserved by the formatter.

## Problem Statement

The linter and language server format Nextflow code by reproducing source text from the AST. A comment can therefore only be printed if it was stored on some AST node during parsing.

Comment capture is opportunistic and incomplete:

- The lexer retypes comments as newline tokens on the default channel. Comments are syntactically newlines, so the parser has no notion of a comment at all.
- The AST builder recovers comments at specific sites (script declarations, statements, etc) by walking the leading/trailing newlines and separators, and stores them as node metadata.
- The formatter emits comments from the node metadata. Any comments not captured in this way cannot be emitted.

Preservation is therefore an allowlist of source positions, not a property of the formatter:

| Input | Current output |
| --- | --- |
| `x = 1  // about x` with a following statement | comment moved to its own line *below* the code |
| `x = 1  // about x` as the last statement | comment deleted |
| comments after the last statement of a block | deleted |
| comments at end of file | deleted |
| comment in an otherwise-empty block | deleted |
| comments inside a list or map literal | deleted |
| comments inside multi-line call arguments | deleted |
| trailing comment on a process output | deleted |

For a tool that rewrites the user's file in place, silent data loss is the worst available failure mode, and it is the most frequently reported formatter problem.

## Goals

- **Formatting never loses, duplicates, or alters the text of a comment.**
- Comments stay in the position a reader would expect: an end-of-line comment stays at the end of its line; an own-line comment stays above the construct it documents.
- Applies to both scripts and config files, and to both formatter front-ends (`nextflow lint -format` and the language server).

## Non-goals

- Preserving original whitespace or indentation; the formatter re-derives indentation.
- New formatting features such as blank-line normalization, length-based line wrapping, or formatting directives. These features may be explored in future efforts.
- Changing Groovydoc behavior. Groovydoc comments (`/** … */`) are stored separately for language-server hover hints. That stays as it is; for *formatting* purposes, a Groovydoc comment is an ordinary leading comment.

## Considered Options

### 1. Total attach pass with a dangling sink

Keep the AST printer; derive all comment attachments in one pass over the comment tokens.

- Good, because the existing printer, options, and tests are reused unchanged -- the diff is confined to comment collection and comment emission.
- Good, because it needs neither a grammar change nor a second lexer pass: comments are already in the parser's token stream, they just weren't all being preserved with the AST.
- Good, because it is the mainstream design. Prettier decorates each comment with `(precedingNode, enclosingNode, followingNode)` and falls back `leading → trailing → dangling(enclosing) → dangling(root)`; Ruff attaches leading/trailing/dangling per node in a side map with a dedicated `place_comment` pass and hand-written overrides.
- Bad, because attachment is a heuristic. Some source positions have no representable emission point, and those comments must be relocated using heuristics.
- Bad, because attachment sites and emission points must be kept in agreement; a comment attached where the printer never asks is a silent loss unless separately asserted.

### 2. Format via `TokenStreamRewriter`

Emit the original token buffer and rewrite only the spans being reformatted.

- Good, because preservation becomes structural. The rewriter makes no modifications to the token stream and only applies queued edits, so any span the formatter does not touch is unchanged.
- Bad, because `WS -> skip` means original indentation is not in the buffer, so every reformatted span still needs the existing pretty-printer.
- Bad, because it requires a much larger change compared to option 1 with no additional benefit.

### 3. On-demand side index (ESLint model)

Attach nothing; let the printer query `commentsBefore/After/Inside` at each emission point.

- Good, because the AST builder's allowlist is no longer needed.
- Bad, because the formatter must query at every emission point. That is the current allowlist problem under a new name.
- Bad, because ESLint is a linter: it reports comment relationships and never has to *decide* where a comment goes. The hard half of the problem is unaddressed.

### 4. Lossless CST (rowan / sprocket model)

Introduce a full-fidelity concrete syntax tree containing every token including trivia, and derive the AST from it.

- Good, because a formatter over a tree that contains every token cannot lose input. Sprocket's `wdl-format` carries comments as first-class trivia (`Trivia::{BlankLine, Comment}`, `Comment::{Preceding, Inline}`) rather than as node attachments, and gets blank-line policy in the same mechanism.
- Bad, because it means maintaining a CST layer in addition to the AST -- a much larger change with new surface area.
- Bad, because `WS -> skip` already forfeits the whitespace half of the fidelity, so the architectural cost buys less than it does for Sprocket.

## Decision

Adopt **option 1**: collect *all* comments from the parser token stream during AST construction, and attach each comment to an AST node during formatting using a set of heuristics.

## Core Capabilities

### Collecting comments from the token buffer

Comments (text and position) are collected from the parser token stream and stored on each module node.

Comments are classified as newlines (`NL`) by the lexer, so they are identified by text prefix -- `//`, `/*`, `#!`. This is unambiguous, because an `NL` token is either a line terminator or a comment. Comments could be given a distinct token type, allowing them to be filtered by type rather than by text prefix, but this change is riskier while providing no additional benefit.

The shebang remains a special case: it is recognized and set explicitly on `ScriptNode`, rather than attached like any other comment. Treating it as an ordinary comment would make it leading on whatever declaration happens to follow it in the source, and declaration sorting could then move it off line 1.

Comment collection is gated on a parser configuration flag: the formatter front-ends enable it, the runtime does not, so a pipeline run does not collect comments for scripts it will never print.

### Attaching comments to nodes

The formatter attaches each comment to an AST node prior to printing. The existing machinery to capture comments during AST construction is removed.

The AST is organized into a tree of **containers**: the module, script declarations with a body, block statements of an `if`/`else`/`try`/`catch`, closures, lists, maps, and argument lists. Each container owns an ordered list of *slots*: an AST node together with the source range in which a comment belongs to it. Containers nest, and sibling containers within a parent are disjoint, so the innermost container holding a comment is found by descent rather than by search.

The reason why a slot can have its own range distinct from the node, is a **method chain**. `a.b().c().d()` is a single nested expression whose sub-expressions all start at `a`, so extents cannot separate one link from the next. Instead, each link gets a slot spanning the gap between the end of its receiver and the start of its method name -- exactly the point where the formatter breaks the line when it wraps the chain. A statement whose root expression is a wrappable chain therefore becomes a container of its own, holding one slot per link. This way, most comments in a method chain are preserved where they are.

Each comment is resolved against its container:

1. **Classify.** *End-of-line* if a non-comment token precedes the comment on the same line; otherwise *own-line*.

2. **Decorate.** Scan the container's slots for `precedingNode` (slot ending before the comment), `followingNode` (first slot starting after it), and `enclosingNode` (innermost slot whose range contains it).

3. **Attach.**
   - end-of-line, and `precedingNode` is on the same line → `TRAILING` on `precedingNode`
   - end-of-line, and `enclosingNode` names a trailing target on the same line (a chain link's receiver) → `TRAILING` on that receiver
   - otherwise `enclosingNode` exists → `LEADING` on `enclosingNode`
   - otherwise `followingNode` exists → `LEADING` on `followingNode`
   - otherwise end-of-line and `precedingNode` exists → `TRAILING` on `precedingNode`
   - otherwise → `DANGLING` on the container, or `LEADING` on it for a container that has no place to print a dangling comment (a chain statement)

Like the allowlist, this approach is fundamentally best-effort. The difference is that it is best-effort in *where comments are placed*, not whether they are placed at all. The attachment rules guarantee that **every comment lands somewhere**. If a comment cannot be preserved exactly where it is -- a comment in the middle of a binary operator expression, for example -- it is hoisted to a leading comment on its enclosing statement (cf. Prettier, Ruff, Black). The user can learn the limitations of the formatter and place comments accordingly, but this is much easier to do when comments are merely reassigned rather than deleted.

This approach also builds on the existing *leading* and *trailing* modes, and adds a third mode -- *dangling* -- which handles comments after the last statement in a container (cf. Ruff).

### Emitting comments

The formatter emits leading, trailing, and dangling comments at every point where a comment can be attached:

- statements and declarations -- leading and trailing
- before a closing brace -- dangling
- end of module -- dangling
- list elements, map entries, call arguments -- leading and trailing
- process and workflow section labels
- method-chain links -- leading and trailing

This logic largely replaces and augments the existing emission logic in the formatter.

The main caveat: **every attachment site must have a matching emission point.** Attachment and emission must be kept in sync in order to guarantee that all comments are preserved. This invariant is verified by testing (see below).

### Verification

The test suite covers all error cases presented in the Problem Statement, as well as all scripts in `docs/snippets` and `tests`. For each test case, the test harness verifies that all comments are preserved and `format(format(x)) == format(x)` (idempotence).

As a final fallback, `nextflow lint -format` refuses skips formatting for a file (and warns) if any comments are missing from the formatter output. This catches alteration and duplication, and it replaces the current silent corruption with a loud, safe failure for every future edge case, including ones not yet imagined.

## Links

- Community issues: nextflow-io/nextflow#6365; nextflow-io/language-server#111, #116, #127, #140
- Prior art:
  - [Ruff](https://github.com/astral-sh/ruff)
  - [Prettier](https://github.com/prettier/prettier)
  - [Sprocket](https://github.com/stjude-rust-labs/sprocket)
  - [ESLint](https://eslint.org/docs/latest/extend/custom-rules)
- Builds on: [Strict syntax parser](20250508-strict-syntax-parser.md)
