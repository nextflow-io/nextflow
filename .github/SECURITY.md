# Security Policy

## Supported Versions

Nextflow uses a calendar-based versioning scheme (`YY.MM.PATCH`).
Releases in April (`04`) and October (`10`) are **stable** releases;
all other months produce **edge** releases (e.g. `26.09.x-edge`).
The currently recommended version is the latest stable release.

| Version stream        | Supported          | Notes                                          |
|-----------------------|--------------------|------------------------------------------------|
| ≤ 24.10.x             | :x:                | End of life                                    |
| 25.04.x               | :x:                | Superseded by 25.10                            |
| 25.10.x               | :white_check_mark: | Previous stable — critical fixes only (\*)     |
| 26.04.x               | :white_check_mark: | Previous stable — critical fixes only (\*)     |
| Edge (`26.xx.x-edge`) | :grey_question:    | Best-effort; may not receive backports         |
| **26.10.x**           | :white_check_mark: | **Currently recommended stable**               |

(\*) Older stable streams are no longer the primary focus of the maintainers,
but critical security fixes are applied when feasible and when reachability
in that stream is confirmed.

## Threat Model

Before reporting, please consult the project [Threat Model](../THREAT_MODEL.md).

Nextflow is a workflow **execution** engine — running the pipeline scripts,
process blocks, configuration closures, and plugin JARs it is given is the
product, not a vulnerability. The threat model states precisely what Nextflow
treats as a security issue and what it does not, so
genuine issues can be distinguished from expected engine behaviour. It also
guides triage of automated code-scanning results, which may arrive as a single
report bundling many candidate findings: each is assigned one disposition from
the threat model, and only confirmed issues enter the private reporting flow
below.

**Compiling untrusted Nextflow DSL is already code execution — not only
running it.** Nextflow's DSL is compiled by `GroovyClassLoader` with full
JVM authority. A service that *compiles* operator-supplied `.nf` scripts
or `nextflow.config` files for analysis, linting, or validation — but
deliberately never executes the resulting classes — is therefore **not** a
safety boundary: code has already run during compilation (static initializers,
`@Grab` dependency resolution, evaluated annotation members, and global AST
transforms on the compile classpath, among other mechanisms). Treat
"compile-only" as equivalent to "run."

See [`THREAT_MODEL.md`](../THREAT_MODEL.md) §§1–3 for the full discussion of
what is in and out of scope.

## Reporting a Vulnerability

Do **not** open a public GitHub issue, pull request, or discussion for a
suspected vulnerability — that defeats coordinated disclosure.

Report it privately using GitHub's built-in
[private vulnerability reporting](https://github.com/nextflow-io/nextflow/security/advisories/new)
feature. This creates a private draft security advisory visible only to
maintainers and the reporter.


When reporting, please include:

- A description of the vulnerability and its potential impact.
- The affected version(s) and configuration.
- Step-by-step reproduction instructions or a minimal proof-of-concept.
- Any suggested mitigations you have identified.

We aim to acknowledge receipt within **3 business days** and to provide an
initial assessment within **10 business days**. Coordinated disclosure is our
default: we will work with you to agree a disclosure timeline before any
public announcement.

## Disclosure Hygiene for Contributors

Until a fix has been publicly announced, do **not** reveal the security nature
of a change anywhere public — commit messages, pull request titles or bodies,
GitHub issues, or review comments — even when the fix touches
security-adjacent code (parameter escaping, deserialization, classloading,
configuration parsing). Describe the behaviour change neutrally. A public
commit or pull request that advertises "fixes the CVE", "security fix", or
"patches the vulnerability" discloses the issue before it is announced and
defeats the coordinated-disclosure process above.

This applies to every contributor, and identically to any AI tooling acting
on a contributor's behalf.
