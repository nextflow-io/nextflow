# no-ext-opts: ADR #7260 vs input-based tool args

Evaluates [nextflow-io/nextflow#7260](https://github.com/nextflow-io/nextflow/pull/7260) (module params + `ToolOpt`/`ToolArgs`, superseding #6769 and #6770) against the input-based approach of [nf-core/methylseq#625](https://github.com/nf-core/methylseq/pull/625) / `../no-ext`, and prototypes the part of the ADR that does help.

Copy of `../no-ext` with one change: tool defaults are **option maps** rendered by `cli()`, and `--args.<tool>` (whole-string replace) became `--opts.<tool>.<option>` (per-option merge). Modules are unchanged (`args: String?`).

```bash
./test.sh                                # 4 scenarios, prints PASS
cd methylseq && nextflow run check.nf    # methylseq defaults: maps vs strings, all param combinations (check.nf itself is untyped test code and fails lint)
```

## ADR recap

- Module-level `params {}` declares tool options (`bwa: ToolArgs = record(K: ToolOpt, Y: ToolOpt)`); `${params.bwa}` renders `-K 100000000 -Y` (null/false omitted, true → bare flag, 1-char → `-`, else `--`).
- Set per included process name: `withName: 'BWA_MEM' { params.bwa.K = ... }` or `-process.BWA_MEM.params.bwa.K=...`.
- Params are resolved once at launch; per-task values are explicitly out of scope ("belongs in a process input").

## Would it work for methylseq?

Only partially, and it would not simplify nf-core/methylseq#625. Every `ext` setting in methylseq, by what it depends on:

| Depends on | Examples | ADR |
|---|---|---|
| sample fields (`single_end`, `id`) | trimgalore `--clip_r2`/3' R2, bismark `--minins/--maxins`, methylation extractor overlap/ignore_r2, every `prefix` (`${id}.sorted`) | ✗ out of scope → still needs `args`/`prefix` inputs |
| samplesheet columns | `trimgalore_args`, `bismark_align_args` | ✗ out of scope |
| pipeline params only | genome prep, methyldackel, qualimap, methurator, multiqc title, coverage2cytosine | only as config (`withName { params.tool.x = params.pbat ? ... }`) |
| constants | markduplicates, addorreplacereadgroups, collecthsmetrics, `--low-memory` | ✓ via `withName` config |
| call site | `SAMTOOLS_SORT` vs dedup sort prefix | ✓ via included name |

Consequences:

- **Two mechanisms per module.** The per-sample cases still need nf-core/methylseq#625's inputs, so trimgalore would get options from `params.trimgalore` *and* `args` — the preset cascade can't be split cleanly because `clip_r2` depends on `single_end`:
  ```groovy
  withName: 'TRIMGALORE' {
      params.trimgalore.clip_r1 = params.clip_r1 > 0 ? params.clip_r1 : params.pbat ? 8 : ...
      params.trimgalore.clip_r2 = ...   // needs single_end: not expressible at launch
  }
  ```
- **Logic moves back into config.** Everything the ADR *can* express is set in `withName` blocks — the config nf-core/methylseq#625 deleted. That reverses "code says *what*, config says *how*", and the config sees untyped `params`.
- **Binding by process name doesn't compose.** methylseq includes `SAMTOOLS_SORT` in both the bismark and bwameth subworkflows; a bare name is ambiguous and a fully-qualified one (`NFCORE_METHYLSEQ:METHYLSEQ:...`) breaks when the pipeline is itself included elsewhere. Args set in workflow code travel with the subworkflow.
- **Module authors will choose inputs anyway** (as argued on #7260): the author can't know whether a caller needs per-sample values, and an input covers both cases.
- **The joins stay.** The complexity nf-core/methylseq#625 added (`ch_meta`/`ch_args` joins, per-subworkflow args records, `toolArgs`) exists for per-sample args, which is exactly the part the ADR excludes.

## What does help: the tool-option model, delivered as inputs

The ADR's serialization rules are useful independently of params. Applied to nf-core/methylseq#625 as untyped maps (works today, `workflows/args.nf`):

```groovy
def cli(opts: Map<String,?>) -> String     // Boolean → flag or omitted, null omitted, prefix inference, `-@` literal keys

def toolArgs(tool: String, s: Sample, p: MiniParams, defaults: Map<String,?>) -> String {
    return s.tool_args[tool] ?: cli(defaults + cliOpts(p.opts[tool] ?: [:]))
}
```

`methylseq/opts.nf` rewrites methylseq's defaults this way; `methylseq/check.nf` asserts identical rendered args for 50,176 param/sample combinations (0 mismatches; a deliberate mutation is caught).

| Function | string version | option map |
|---|---|---|
| `trimgaloreArgs` | 61 lines (nested ternaries ×4) | 20 (preset table + one line per option) |
| `bismarkAlignArgs` | 22 | 20 |
| `bismarkMethylationExtractorArgs` | 12 | 14 |

So the code win is local: large where there is a preset cascade, neutral elsewhere (`cond ? '--x' : ''` becomes `x: cond`). The user-facing win is bigger:

- `--opts.bismark_align.maxins=700` changes one option; `--args.bismark_align='...'` had to restate every default.
- `--opts.trimgalore.fastqc=false` removes a default flag.
- Unlisted options pass through (`--opts.trimgalore.quality=30`) — no enumeration needed (pinin4fjords' `**kwargs` concern on #6770).
- Nested CLI maps parse as-is: `opts: Map<String,Map<String,?>>`. CLI values arrive as strings, so `cliOpts()` turns `'true'`/`'false'` into Booleans; string values in defaults (picard's `ASSUME_SORTED: 'true'`) render verbatim.

Ceilings:

- Values are `?` to the type checker (`Map<String,?>`); the linter rejects `Map.findAll`/`collectEntries`/`Entry.key` and `?` as a standalone type, hence `keySet()`, `inject` and `instanceof Boolean`. No per-option typing — that is what a real `ToolOpt`/`ToolArgs` type would add.
- User opts bypass sample-dependent gating (`--opts.align.maxins` also reaches single-end samples; see `test.sh`).
- Not representable: repeated flags (`-I a -I b`), positional args, `--k=v` / glued `-O2` separators. User-added options render after defaults (map insertion order).
- The samplesheet `<tool>_args` column is still a raw string that replaces the whole thing; a raw string can't merge with a map.
- nf-schema validation of a nested `object` param untested (the example doesn't use nf-schema).

## Recommendation

- Keep nf-core/methylseq#625's delivery (inputs set in workflow code). Module params don't cover methylseq's per-sample args, and the rest would move logic back to config.
- Option maps + `cli()` + `--opts` are a cheap, optional follow-up to nf-core/methylseq#625: rewrite `args.nf` defaults, replace `params.args` with `params.opts`, modules untouched. Worth it mainly for per-option overrides. Ported on the `refactor-ext-meta` branch as a follow-up commit; every task script in the 9 verification configs matches the string-based version ignoring whitespace.
- User option values render verbatim: `--opts.multiqc.title='"My run"'` needs its own quotes (the pipeline default quotes `--multiqc_title`).
- Feedback for #7260: the valuable part is the `ToolArgs`/`ToolOpt` *type* (documented options + rendering rules). Making it usable as a process input type — `input: record(id: String, reads: Path, args: BwaArgs?)` with `bwa mem ${args}` — would give typed, documented, per-task options, and `nextflow module run` could set it like any other input. That covers both problems stated on #7260 (documenting options, supplying args to `module run`) without a second params namespace or name-based binding.
