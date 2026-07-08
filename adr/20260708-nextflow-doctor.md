# Diagnose broken Nextflow installs with `nextflow doctor`

- Authors: Edmund Miller, OpenAI
- Status: draft
- Date: 2026-07-08
- Tags: cli, install, diagnostics, troubleshooting

## Summary

Add a `nextflow doctor` command with a Bash bootstrap fallback plus a core JVM command. The fallback keeps diagnostics available when Java is missing, Java is unsupported, or the framework download is broken; the JVM command handles install checks once Nextflow can start normally.

## Problem Statement

The `nextflow` launcher checks Java and downloads the framework before Groovy commands run. When either bootstrap step fails, users get a one-off fatal message instead of a consolidated install health report.

That leaves users debugging install problems by trial and error before they can even run `nextflow info` or a pipeline. Nextflow should provide one diagnostic entry point that covers both pre-JVM bootstrap failures and JVM-visible install problems.

## Goals

- Diagnose Java availability/version before the JVM command starts.
- Diagnose release server/framework download reachability.
- Diagnose local `$NXF_HOME`, temp, config, registry, and Seqera Platform endpoint health once the JVM starts.
- Produce stable text output and non-zero exit status when required checks fail.
- Avoid downloads, repairs, credentials disclosure, and pipeline execution.

## Non-goals

- No repair mode.
- No JSON output in this first cut.
- No container, executor, or cloud-provider validation.
- No project-specific pipeline validation.

## Decision

### D1 — Hybrid command surface

Add `nextflow doctor` as both a Bash launcher fallback and a core JVM command.

The Bash launcher handles bootstrap-only diagnostics when Java is missing, Java is unsupported, or the framework JAR is absent. This path runs before the JVM command can start and must stay small.

The core `CmdDoctor` command handles normal diagnostics after Java and the framework are available.

### D2 — Stable text report and exit status

The default report is stable, plain text with uncolored status lines:

```text
Nextflow doctor
[OK] Java runtime - Java 21.0.4
[FAIL] Release server - https://www.nextflow.io/releases/latest/version?current=25.04.0
Doctor summary: 1 OK, 0 WARN, 1 FAIL, 0 SKIP
```

`FAIL` means the command exits non-zero. `WARN` and `SKIP` do not make the command fail.

### D3 — Diagnostic-only bootstrap

The bootstrap fallback must not download the framework JAR, repair files, or write Java version cache entries.

`nextflow doctor` reports install state. It does not mutate the install.

### D4 — Install-scoped JVM checks

The JVM command checks:

- Java runtime;
- Nextflow runtime;
- Nextflow home directory;
- temporary directory;
- default configuration parsing;
- release server reachability;
- plugin registry reachability;
- Seqera Platform reachability.

It does not validate containers, executors, credentials, cloud accounts, SCM providers, or pipeline scripts.

### D5 — Fast endpoint probes

Network checks use one request with a fixed 5 second timeout.

The doctor command is a health check, not a retrying updater. Slow retry policies hide the failing endpoint and make troubleshooting worse.

### D6 — Offline behavior

When `NXF_OFFLINE=true`, network probes are skipped and reported as `SKIP`.

Offline mode should still report local Java, runtime, home, temp, and configuration health.
