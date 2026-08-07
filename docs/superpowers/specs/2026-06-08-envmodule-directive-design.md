# Spec: Rename `module` process directive to `envModule`

- **Date:** 2026-06-08
- **Status:** Draft
- **Author:** Paolo Di Tommaso

## Summary

Rename the Environment Modules process directive from `module` to `envModule` to
disambiguate it from the Nextflow *module* concept (script `include` and the
module registry system under `nextflow/module/*`). The old `module` directive is
**retained as a deprecated alias** so existing pipelines keep working; using it
emits a runtime warning advising migration to `envModule`.

This is a non-breaking change: no behavior changes for users, only a new
preferred name and a deprecation warning for the old one.

## Motivation

The `module` directive configures [Environment Modules](http://modules.sourceforge.net/)
required by a task. The name collides conceptually with two unrelated Nextflow
"module" features:

- **Script modules** — reusable libraries pulled in via `include { ... } from`.
- **The module registry** — the `nextflow module` CLI and the `nextflow/module/*`
  packages for installing/resolving versioned modules.

The overloaded term confuses users and makes documentation harder to write.
Renaming the directive to `envModule` makes its purpose explicit while preserving
backward compatibility.

## Goals

- Introduce `envModule` as the preferred directive for Environment Modules.
- Keep `module` working as a deprecated alias.
- Emit a runtime deprecation warning when `module` is used, pointing to `envModule`.
- Update reference documentation to lead with `envModule`.

## Non-goals

- No changes to the Nextflow module/registry system or `include`.
- No removal of the `module` directive (deprecation only — removal is a future,
  separate decision).
- No changes to how Environment Modules are loaded at runtime.

## Design

### Reference precedent

The `echo` → `debug` rename is the established pattern in this codebase. `echo`
remains as a method that logs a deprecation warning and forwards to `debug`:

```groovy
// ProcessBuilder.groovy
void echo( value ) {
    log.warn1('The `echo` directive has been deprecated - use `debug` instead')
    config.put('debug', value)
}
```

`envModule` follows the same shape.

### Internal config key: keep `module`

Both `envModule` and the deprecated `module` resolve to the **same internal
config key**, which remains `module`. This is deliberate:

- `TaskConfig.getModule()`, `TaskBean.moduleNames`, `TaskHasher`, and
  `BashWrapperBuilder` continue to read the `module` key unchanged.
- Task hashing is unaffected, so this rename does **not** invalidate existing
  resume caches.
- The change is confined to the DSL surface (the directive method names and the
  whitelist), keeping the blast radius minimal.

### Changes

1. **`modules/nextflow/src/main/groovy/nextflow/script/dsl/ProcessBuilder.groovy`**
   - Add `void envModule(value)` containing the current `module()` implementation
     (append to the `ConfigList` stored under the `module` key).
   - Replace the body of `void module(value)` to emit
     `log.warn1('The `module` directive has been deprecated - use `envModule` instead')`
     and delegate to `envModule(value)`.
   - Add `'envModule'` to the directives whitelist. Keep `'module'` in the list so
     the alias remains valid.

2. **`modules/nextflow/src/main/groovy/nextflow/processor/TaskConfig.groovy`**
   - The `put()` special-case for dynamic list values currently keys off
     `'module'`. Since the internal key stays `module`, no change is required
     here. (Confirm during implementation that `envModule` values still land
     under the `module` key before this check runs.)

3. **`modules/nf-lang/src/main/java/nextflow/script/dsl/ProcessDsl.java`**
   - Add an `envModule(String value)` declaration with the descriptive doc
     comment currently attached to `module`, pointing at the updated docs anchor.
   - Keep the existing `module(String value)` declaration so the language server
     still recognizes the alias. (Per decision, no `@Deprecated` annotation /
     language-server diagnostic in this change — runtime warning only.)

4. **`docs/reference/process.md`**
   - Rename the `### module` section to `### envModule` (update the
     `(process-module)=` anchor accordingly, e.g. `(process-envmodule)=`, and add
     a redirect/alias anchor if needed to avoid breaking inbound links).
   - Update examples to use `envModule`.
   - Add a short note: "The `module` directive is a deprecated alias for
     `envModule` and will emit a warning when used."

### Warning behavior

- The warning is emitted via `log.warn1` (de-duplicated — logged once per run)
  at process-build time when `module` appears in a process definition or config.
- Message: `The `module` directive has been deprecated - use `envModule` instead`.

## Backward compatibility

- Existing pipelines using `module` continue to run unchanged.
- Resume caches remain valid (internal key unchanged → identical task hashes).
- Config files (`nextflow.config`) setting `process.module` continue to work; the
  deprecation surfaces through the same DSL path.

## Testing

- **Unit (Spock):**
  - `envModule` populates the same config as `module` (single value, list, and
    `:`-separated forms via `TaskConfig.getModule()`).
  - `module` still works and triggers the deprecation warning (assert via log
    capture, mirroring existing `echo`/`debug` tests if present).
  - `envModule` is accepted as a valid directive; unknown directives still rejected.
- **Hashing:** a task using `envModule 'x'` produces the same hash as the same
  task using `module 'x'` (guards against cache invalidation).
- **Docs snippet:** ensure the updated `envModule` example in `process.md` is valid.

## Rollout

- Ship in an edge release first.
- Document in `CHANGELOG.md` under the next version.
- Removal of the `module` alias is explicitly deferred and not scheduled by this spec.
