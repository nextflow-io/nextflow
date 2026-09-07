# Agents as modules

An `agent` can be authored as its own **module** and consumed with the ordinary
`include` statement, exactly like a process — including under an alias. The module
directory carries the agent's own skills and its own tool, so it is a shippable
unit rather than a snippet to copy into a script.

## What's here

```
main.nf                                   include { reporter as qc } from './mods/reporter'
mods/reporter/
  main.nf                                 the agent: skills 'qa-report','style' + tools 'nf:module_run'
  skills/qa-report/SKILL.md               report format (module-local skill)
  skills/style/SKILL.md                   house style   (module-local skill)
  tools/qc_verdict.nf                     the deterministic PASS/WARN/FAIL rule, included by the module
```

Three things are being demonstrated at once:

1. **The include, aliased.** `include { reporter as qc } from './mods/reporter'`
   resolves the directory to `mods/reporter/main.nf`. The alias renames the task,
   so the work dir, progress table and trace all say `qc`.
2. **Module-local skills.** Both skills resolve under
   `mods/reporter/skills/<name>/` — the directory of the file that *declares* the
   agent, **not** the directory of the script that includes it, and that holds
   under the alias. There is no fallback to the project dir, so a
   `skills/qa-report/` next to `main.nf` could not shadow the module's own.
3. **Module-scoped tools.** `tools 'nf:module_run'` sees only what
   `mods/reporter/main.nf` itself defines or includes — hence the module includes
   its own `tools/qc_verdict.nf`. A process defined by the *including* script is
   deliberately invisible, so the agent's tool surface does not change depending
   on who imported it.

The `qc_verdict` tool also splits the work honestly: the verdict is computed by
code — a plain `script:` block whose `stdout()` becomes the tool result — and the
prose is the model's job.

## Run it

```bash
export OPENAI_API_KEY="sk-..."
nextflow run main.nf
```

Expected:

```
ANSWER=
[MOD-QA v1]
VERDICT: PASS
METRICS: N50 = 45 kb, completeness = 96.4%, total length = 5.1 Mb
NOTE: Metrics clear the thresholds; proceed. -- mods/reporter
```

`[MOD-QA v1]` and the `-- mods/reporter` sign-off are the tell that both
module-local skills were activated; `VERDICT` comes from the tool, not the model.
Add `-with-agent-trace` to watch the `activate_skill` and `qc_verdict` calls.

To prove the skills really come from the module directory, create
`skills/qa-report/SKILL.md` next to this `main.nf` with different instructions and
re-run: the output does not change.

## Notes

- **`nextflow.enable.types` is per file.** An `agent` block does not need it; the
  typed process in `tools/qc_verdict.nf` does. Neither this script nor the agent
  module sets it.
- **Config selectors work from the outside** (see `nextflow.config`): an
  `agent { withName: … }` block matches the declared name `reporter`, the alias
  `qc`, or the fully-qualified name — so the consumer can place the agent without
  editing the module. Agents read the `agent` scope only; a `process` selector
  matching an agent's name has no effect on it.
- **Params are inherited** from the run; a module agent's prompt may read
  `params.*`. `addParams`/`params` on the include are not an agent feature.
- **Resume**: editing a module skill invalidates the agent's cache entry; moving
  or renaming the module directory does not; aliasing does (the task name is part
  of the key).
- **Registry-hosted agent modules** (`include { a } from 'scope/name'`) are not
  supported yet — local paths only.
