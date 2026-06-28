# Agent skills

Demonstrates the `skills` directive: giving an agent one or more Anthropic-style
**skills** (`SKILL.md` folders) that progressively disclose instructions to the
model. Skills extend an agent's capability the way `tools` do, but instead of
calling code they inject expert instructions on demand.

## What's here

```
skills/
  main.nf                         # an agent declaring `skills 'sequence-report'`
  skills/sequence-report/SKILL.md # a local skill: a standardized QC report format
```

The `reporter` agent declares `skills 'sequence-report'`. That name resolves to
the local `skills/sequence-report/` directory (the `skills/` directory sits
beside `main.nf`). At run time the model sees the skill's `name` + `description`
in an available-skills catalog; when the prompt asks for an assembly summary it
calls `activate_skill` to read the full instructions, then formats its answer as
the skill dictates.

## Run it

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected: the answer is the distinctive skill-dictated format, e.g.

```
ANSWER=
[SEQ-REPORT v1]
STATUS: PASS
METRICS: N50 = 45 kb, total length = 5.1 Mb, GC = 50.8%
NOTE: A 5.1 Mb assembly with 45 kb N50 looks suitable to proceed.
```

The `[SEQ-REPORT v1]` header is the tell that the skill was activated — without
the skill the model would answer in free prose. Add `-with-agent-trace` to see
the `activate_skill` call in the execution trace.

## Local vs remote skills

- **Local** — a bare name (e.g. `sequence-report`) resolves to
  `skills/<name>/` beside the script.
- **Remote** — a GitHub reference (`github.com/<org>/<repo>[@rev]`,
  `https://github.com/...`, or `git@github.com:...`) is cloned and cached into
  the same `skills/` directory, then loaded the same way.

> **Trust note:** a remote skill's `SKILL.md` becomes model instructions, so
> activating a remote skill means trusting its authors. Pin a commit SHA
> (`@<sha>`) rather than a moving branch so the content can't change under you.
