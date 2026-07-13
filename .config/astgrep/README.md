# ast-grep lint rules for Groovy

Structural (AST-based) lint rules for Nextflow's Groovy sources, using
[ast-grep](https://ast-grep.github.io/) with a compiled tree-sitter Groovy grammar
registered as a custom language.

## Usage

```bash
cd .config/astgrep
./build.sh                                    # fetch pinned grammar + compile groovy.dylib (run once)
cd -
ast-grep scan -c .config/astgrep/sgconfig.yml # scan the repo
```

Requires [`ast-grep`](https://ast-grep.github.io/guide/quick-start.html) (>= 0.44) and a
C compiler (`gcc`/`clang`).

## How it works

- **`build.sh`** fetches [dekobon/tree-sitter-groovy](https://github.com/dekobon/tree-sitter-groovy)
  at a pinned immutable commit (`8c70dc6`, release v0.2.2) into `.grammar-build/` (gitignored,
  not vendored), asserts the resolved SHA matches the pin, and compiles the parser into
  `groovy.dylib`. The library is native code compiled on each machine, so it works on macOS
  (Mach-O) and Linux (ELF) alike — no per-platform binary is committed.
- **`sgconfig.yml`** registers `groovy` as an ast-grep custom language (`expandoChar: _`
  because Groovy uses `$` for GString interpolation; metavariables are still `$NAME` / `$$$`).
- **`rules/`** holds one YAML lint rule per file.

`groovy.dylib` and `.grammar-build/` are gitignored; each checkout rebuilds the library locally.

## Rules

| id | severity | flags |
|---|---|---|
| `no-wildcard-import` | warning | `import foo.bar.*` — import specific classes |
