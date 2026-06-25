# Separate CLI from Runtime

- Authors: Ben Sherman
- Status: accepted
- Deciders: Ben Sherman, Paolo Di Tommaso
- Date: 2026-06-25
- Tags: cli, modularization, build, architecture

## Summary

Move the CLI code into a new `nf-cli-v1` module, leaving `nextflow` as a pure runtime that knows nothing about the CLI. The split untangles circular dependencies, speeds up builds, lets library consumers depend on the runtime without dragging in the CLI, and paves the way for a future CLI v2.

## Problem Statement

The `nextflow` module mixed two distinct concerns:

1. **The CLI**: the `Launcher` entry point, argument parsing (JCommander), and the `Cmd*` command classes (`CmdRun`, `CmdConfig`, `CmdLog`, the `module` subcommands, etc.).
2. **The runtime**: the dataflow engine, session, executors, config model, SCM/asset management, and everything required to actually execute a workflow.

Bundling them together caused several problems:

- **Circular dependencies**: CLI code reached into the runtime and runtime code reached back into CLI concerns (e.g. `ConfigBuilder` carried CLI-specific logic), making the dependency graph hard to reason about.
- **Slow builds**: editing a CLI class forced a rebuild of the entire module, including the runtime.
- **Unwanted coupling for consumers**: library consumers -- Seqera Platform, plugins -- that only need the runtime were forced to depend on the whole CLI surface.
- **No room for a CLI v2**: a redesigned CLI could not be introduced cleanly while the existing CLI was fused to the runtime.

In addition, several plugins (`nf-console`, `nf-k8s`, `nf-tower`, `nf-wave`) contribute CLI commands, so naively moving the CLI out risked coupling those plugins to the new CLI module instead of to the runtime.

## Goals or Decision Drivers

- **Separation of concerns**: the CLI should depend on the runtime, never the reverse.
- **Faster, more incremental builds**: changing CLI code should not rebuild the runtime.
- **Runtime reusable as a library**: consumers (Platform, plugins) can depend on `nextflow` alone, without the CLI.
- **Break circular dependencies**: produce a clean, acyclic module graph.
- **Keep plugins decoupled from the CLI**: plugins that add commands should depend on the runtime, not on `nf-cli-v1`, wherever feasible.
- **Enable a future CLI v2**: make it possible to add an alternative CLI implementation alongside v1.
- **Preserve user-facing behavior**: the Nextflow CLI should continue to behave the same way.

## Non-goals

- **Redesigning the CLI**: this is a structural move only; CLI v2 is explicitly future work. The module is named `nf-cli-v1` to leave that door open.
- **Removing the `nf-cli-v1` → `nf-k8s` build-time dependency**: `K8sDriverLauncher` now lives in `nf-cli-v1` and pulls in `nf-k8s` at build time (see below). Loading `nf-k8s` on demand instead -- so the CLI need not depend on it at build time -- is left to CLI v2.
- **Changing the packaging/distribution format**: the produced `nextflow` binary is unchanged from the user's perspective.

## Solution

Extract the CLI into a new `nf-cli-v1` module that depends on the `nextflow` runtime, and relocate the shared interfaces and adapter classes needed to keep the dependency graph acyclic and plugins decoupled from the CLI.

## Rationale & discussion

### New module structure

```
nf-cli-v1 ──► nextflow (runtime) ──► nf-commons, nf-httpfs, nf-lang
   │
   ├──► nf-k8s (via CmdKubeRun) ──► nextflow
   └──► nf-lineage (via CmdLineage) ──► nextflow
```

- `nf-cli-v1` contains the `Launcher` entry point, `CliOptions`/`HubAware`, and all `Cmd*` classes (including the `module` subcommands). Its `application` main class is `nextflow.cli.Launcher`, and it produces the shadow jar that becomes the `nextflow` distribution.
- `nextflow` contains the runtime only -- no CLI entry point, no `Cmd*` classes.
- `settings.gradle` registers the new `nf-cli-v1` module.
- `packing.gradle` now packs from `:nf-cli-v1:shadowJar` instead of `:nextflow:shadowJar`; the produced `nextflow-<version>-one.jar` / `-dist` artifacts are unchanged in name and behavior.

### Breaking the circular dependencies

Several deliberate moves were required so that the runtime never depends on the CLI:

1. **`ConfigBuilder` → `ConfigCmdAdapter`**: the CLI-specific portions of `ConfigBuilder` were extracted into a new `ConfigCmdAdapter` class in `nf-cli-v1`, leaving the core config builder in the runtime. This also yields a cleaner separation of concerns.

2. **Command interfaces → runtime**: `AuthCommand` and `LaunchCommand` interfaces were moved into the runtime. This lets `nf-tower` implement these commands while depending only on the runtime and not the CLI. `nf-lineage` was similarly reworked to depend on `nextflow` and be required by `nf-cli-v1`, removing the need for the `LinCommand` interface (the lineage command uses `LinCommandImpl` directly).

3. **`PluginExecAware` → runtime**: moved into the runtime (`nextflow.plugin.cli`) so that plugins can declare CLI commands without depending on `nf-cli-v1`. The `exec()` method was changed in order to remove the dependency on `nextflow.cli.Launcher`.

### Plugins that contribute CLI commands

`nf-console`, `nf-k8s`, `nf-tower`, and `nf-wave` add CLI commands. The goal was to keep them depending on the runtime rather than the CLI:

- `nf-console` still depends on `:nf-cli-v1`.

- `nf-k8s` compiles against `nextflow` only and references no CLI classes. To get there, `K8sDriverLauncher` (which uses `CmdRun`) was moved out of `nf-k8s` and into `nf-cli-v1`. As a result, `nf-cli-v1` now depends on `nf-k8s` at build-time. The CLI v2 will not need to do this if it does not preserve the `kuberun` command; it can simply load `nf-k8s` at runtime when using the `k8s` executor.

- `nf-tower` and `nf-wave` were decoupled by moving the relevant command interfaces into the runtime; both now compile against `nextflow` only.

### Breaking change: plugin commands

The `exec()` method of `PluginExecAware` was changed to remove the dependency on the CLI:

```groovy
// before
int exec(Launcher launcher, String pluginId, String cmd, List<String> args)

// after
int exec(String pluginId, String cmd, List<String> args)
```

As a result, `PluginAbstractExec` (the base class used by most plugin commands) no longer loads the Nextflow configuration the way it used to -- previously it ran `ConfigBuilder` over the launcher options and config files. It now creates the `Session` from a **minimal config** derived from environment variables only:

- `workDir` from `NXF_WORK` (default `work`)
- `cloudcache` from `NXF_CLOUDCACHE_PATH` (when set)

There are two plugin commands among the core plugins that are affected:

- **`CacheCommand`** (`nf-tower`, the `cache-backup`/`cache-restore` commands) -- adapted to read `NXF_CLOUDCACHE_PATH` from the environment directly. The `cloudcache` entry in the minimal config exists specifically for this purpose. The `cache-backup` command is used by Seqera Platform to upload logs to the cloudcache on exit.

- **`WaveCmdEntry`** (`nf-wave`) -- not adapted. These commands construct a `WaveClient` from `session.config.wave`, `session.config.fusion`, and `session.config.tower`, which the minimal config no longer populates. As a result they now honor only environment variables and defaults, **not** settings from `nextflow.config`. These commands are used only for debugging, so they are not critical, but may warrant further review.

Plugin commands that need configuration should load it explicitly rather than relying on `PluginAbstractExec`. This change is acceptable because plugin commands are quite rare in practice.
