# Package provider plugins

- Authors: Evan Rees
- Status: proposed
- Deciders: Nextflow maintainers
- Date: 2026-07-07
- Tags: plugins, packages

## Summary

Provide a single generic `package` process directive whose package managers are implemented as pluggable providers, each shipped as its own `nf-<provider>` plugin.

## Problem Statement

Nextflow supports software dependencies through per-manager process directives (`conda`, `spack`), and each new package manager requires its own directive and core runtime changes. Requests for additional managers keep arriving (nix and guix in #1526, uv in #7023, pixi and R managers in the #6342 review), and the per-manager approach does not scale: every addition grows the core, the process DSL, and the config surface.

## Goals or Decision Drivers

- One directive for all package managers, so processes declare dependencies the same way regardless of the tool.
- New managers must be addable without core runtime changes.
- Alignment with Wave: providers that Wave can build (conda-based) hand off to Wave; others build locally.
- The feature must be opt-in (preview flag) and must not affect existing pipelines.

## Non-goals

- Implementing Wave container builds for non-conda providers.
- Removing or deprecating the legacy `conda` and `spack` directives.
- Extracting the existing conda support out of core (the legacy `conda` directive still depends on it).

## Considered Options

### Per-manager directives

Continue adding one process directive per package manager (`uv`, `pixi`, `nix`, ...).

- **Pro:** Explicit and familiar; mirrors the existing `conda`/`spack` directives.

- **Con:** Every manager requires parser, DSL, config, and runtime changes in core; N managers means N directives to document and maintain.

### Wave as the single integration point

Route all package specifications through Wave and let it build the environments.

- **Pro:** One container-based path for every manager.

- **Con:** Wave only builds conda-based package environments; forcing other managers through it is not possible today, and it makes containers mandatory where local environments suffice.

### Generic directive with pluggable providers (chosen)

A `package` directive backed by a `PackageProvider` SPI in core, with each manager implemented in its own `nf-<provider>` plugin.

- **Pro:** Managers ship, version, and test independently; third parties can add providers without touching core.

- **Con:** `package` is a Groovy/Java keyword, so the strict (v2) parser needs a narrowly scoped exception, and the directive cannot be offered under the legacy parser.

- **Con:** A new SPI surface in core to keep stable.

## Solution

Add a provider-agnostic SPI to core (`PackageProvider`, `PackageProviderExtension`, `PackageManager`, `PackageSpec`, and the `packages` config scope) and implement each manager as a plugin: `nf-conda`, `nf-pixi`, `nf-uv`, `nf-nix`, `nf-guix`, `nf-pak`, `nf-install2r`.

- The directive is gated behind the `nextflow.preview.package` feature flag and requires the v2 syntax parser; the `package` keyword is allowed only in the process-directive section, and the v1 parser is untouched.

- When Wave is enabled, conda-family and pixi specifications map to Wave's CONDA build type; local-only providers skip Wave and their plugin creates the environment on the local file system, so mixed pipelines work.

- When a process has no `package` directive, the process module directory is scanned for a provider manifest file (e.g. `environment.yml`, `requirements.txt`), analogous to Wave's Dockerfile auto-detection; this can be disabled with `packages.autoDetect = false`.

- Conda's environment cache remains in core because the legacy `conda` directive depends on it; the `nf-conda` plugin reuses those classes. Full extraction into the plugin is a follow-up once the legacy directive is deprecated.
