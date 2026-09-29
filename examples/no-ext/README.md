# no-ext: tool args and metadata without `ext` or meta maps

Minimal pipeline exercising the core nf-core/methylseq patterns (trim → align → sort → optional dedup → multiqc) with static typing and records, but **no `ext` config and no `meta` map**. Processes are stubs that echo their command line, so the resolved args can be checked.

Requires Nextflow 26.09.1-edge or later.

## Usage

```bash
nextflow lint .
nextflow run . --input data/samplesheet.csv --fasta data/genome.fa -output-dir results
./test.sh    # runs 3 scenarios and asserts the resolved args
```

Overriding args:

```bash
# per run, for all samples (use '=' when the value starts with --)
nextflow run . ... --args.trimgalore='--quality 30' --args.samtools_sort='-m 2G'
```

```csv
# per sample, via a <tool>_args column
sample,fastq_1,fastq_2,clip_r1,trimgalore_args
s3,data/s3_1.fq.gz,data/s3_2.fq.gz,,--custom-trim
```

Precedence: `<tool>_args` samplesheet column > `--args.<tool>` > pipeline default.

## Layout

| File | Role | Replaces |
|------|------|----------|
| `main.nf` | params, samplesheet → flat records, publish/output | `params` in config, `publishDir` |
| `workflows/mini.nf` | workflow logic; builds `args`/`prefix` per call | `ext.args` / `ext.prefix`, `withName` selectors |
| `workflows/args.nf` | typed functions computing default args from `(sample, params)` | `conf/modules/*.config` closures |
| `workflows/types.nf` | `Sample`, `MiniParams` record types | meta map |
| `modules/*.nf` | processes with `args: String?`, `prefix: String?` record inputs | `task.ext.args`, `task.ext.prefix` |

## Patterns

- **Samplesheet columns are record fields.** `record(id, single_end, reads, clip_r1, tool_args)`. No `meta`, so no need for both `id` and `meta.id` — joins use `by: 'id'`.
- **ext settings are process inputs.** Modules declare `args: String?` and `prefix: String?` in their input record. Nullable fields may be absent from the incoming record.
- **Defaults are typed functions.** `trimgaloreArgs(s: Sample, p: MiniParams)` has the same inputs the config closure had (`meta` + `params`) but is type-checked.
- **One override param, not one per tool.** `params.args: Map<String,String>` is set from the CLI with dotted keys (`--args.<tool>`). Per-sample overrides come from `<tool>_args` columns, collected into `Sample.tool_args` at parse time. `toolArgs(tool, s, p, default)` resolves them.
- **Samplesheet columns falling back to params** are resolved once at parse time: `clip_r1: row.clip_r1 ? row.clip_r1.toInteger() : params.clip_r1`. Downstream code only sees `s.clip_r1`.
- **Per-call-site settings** (formerly `withName`) are just different arguments at each call, e.g. `SAMTOOLS_SORT` vs `SAMTOOLS_SORT_DEDUP` prefixes.
- **Recovering sample fields after a process**: process outputs carry `id` + new files; `ch_samples.join(ch_out, by: 'id')` brings back `single_end` etc.

## Findings

- **Forwarding the input record is not type-safe.** With `input: sample: Sample` and `output: sample + record(...)`, extra fields survive at runtime, but the output type is only the declared fields, so `r.extra` downstream is a lint error. Joining by `id` keeps full types.
- **Stale `args` can leak.** Every module uses the same field names (`args`, `prefix`). If a record carrying `args` is passed to another process, that process silently picks it up. Set `args` only in the `.map` that feeds a call, and never emit it from a module. (Inferred from record merge semantics, not tested.)
- **CLI values starting with `--` need `=`**: `--args.align='--foo'`, otherwise `--foo` is parsed as a param. Existing Nextflow behavior.
- **What users lose vs `ext`:** only tools the developer routes through `toolArgs` are overridable; keys are developer-chosen names, not `withName` glob patterns. A generic "any param overridable per sample" would need a map again — keep per-sample columns explicit.
- **Non-data settings can stay in config** (e.g. `ext.singularity_pull_docker_container`). A GPU flag is a Boolean passed from the entry workflow.
- `collect()` yields `Bag<Path>`, so process inputs must be `Bag<Path>`, not `Set<Path>`.

## Port to nf-core/methylseq

[nf-core/methylseq#625](https://github.com/nf-core/methylseq/pull/625) (branch `refactor-ext-meta`, stacked on `preview-26-04`) removes `ext` config and meta maps from the whole pipeline.

### What changed

- **Modules (40)**: `meta: Record` → `id: String` (+ `single_end: Boolean`, `strandedness: String?` where used) and `args: String?` / `prefix: String?` / `suffix: String?` / `args2: String?` record fields; outputs drop `meta`. Scripts reassign the inputs (`args = args ?: ''`, `prefix = prefix ?: "${id}"`) so the diff per module is small. Modules without a sample record (genome prep, faidx, untar, gunzip, bwa index, sequence dictionary, interval list) take a positional `args: String?` input; their never-configured `ext.prefix` became a fixed default.
- **Config**: all of `conf/modules/` and `conf/subworkflows/` deleted. The one non-ext setting (`cache = 'lenient'` for `PARABRICKS_FQ2BAMMETH`) moved to `conf/base.config`; unused `ext.use_gpu` removed. `ext.singularity_pull_docker_container` (container directive) is left alone.
- **Defaults**: `workflows/methylseq/args.nf` holds every former `ext.args` default as typed functions, plus small params record types (`BismarkParams`, `MethyldackelParams`, ...) documenting which params each default reads.
- **Subworkflow args API**: each nf-core subworkflow declares a per-sample args record in its input type — `bismark_args: BismarkArgs?`, `bwameth_args`, `bwamem_args`, `methyldackel_args`, `targeted_args` — and maps each field to `args` only when calling that step. Run-level tools (genome prep, faidx, dict, intervallist) take plain `String` subworkflow inputs. Structural prefixes (`.sorted`, `.deduplicated.sorted`, `.markdup.sorted`, `targeted.cov`) are set inside the subworkflows.
- **Records**: `Sample = {id, single_end, reads, tool_args}`; `Alignment` / `AlignedSample` are just `{id, bam[, bai]}`. Subworkflows keep a `ch_meta` (`id`, `single_end`) and `ch_args` (`id`, `<sub>_args`) and join them by `id` where needed; `METHYLSEQ` does the same with `ch_meta` for post-alignment tools.
- **User overrides**: `--args.<tool>` for 25 tools (keys listed in the `args` help text in `nextflow_schema.json`), and `trimgalore_args` / `bismark_align_args` samplesheet columns.

Size vs `preview-26-04`: modules +342/−299, conf −339, subworkflows +211/−69, workflows +379/−32 (`args.nf` is 274 lines, ~70 of them params record types).

### Verification

`-profile test,docker`, 9 param sets covering bismark (default; pbat/nextseq/length/mismatches/minins/unmapped/3' clip; zymo+hisat+local+maxins; cytosine_report+skip_dedup+targeted+comprehensive; methurator+qualimap+preseq+title; combined_index), bwameth (targeted+hsmetrics+qualimap+methyldackel opts; skip_dedup+merge_context+ignore_flags), and bwamem+taps. Compared the `.command.sh` of every task and the published file tree against `preview-26-04`:

- Output trees identical for all 9; same tasks succeed/fail (preseq and one bwameth methyldackel config fail on the tiny test data in both).
- All task scripts byte-identical except: MULTIQC (dropped the always-empty `${prefix}` line), `bismark2summary` BAM order (Set order), rastair tasks (untagged, so `(1)`..`(4)` map to different samples; identical as a set), and which tasks ran before the failing config aborted.
- Override run: `--args.fastqc`, `--args.bismark_methylationextractor`, `--args.multiqc` all reach the tool.
- 4 `-stub` runs (bismark, bwameth+targeted, bwamem+taps, qualimap/preseq/methurator) pass.

### Findings

- **Join back metadata only, not the whole input record.** Joining the full reads record after a process dragged `reads` and args fields into every downstream result. Join `record(id, single_end)`; join the args record only at the call sites that need it.
- **Namespace per-sample args records** (`bismark_args`, not `args`): a field named `args` or `prefix` on a pass-through record is silently picked up by the next process that declares it. Attach `args`/`prefix` only in the `.map` feeding a call; modules never emit them.
- **One nested args record per subworkflow scales better than one field per step**: callers build it with one function (`bismarkArgs(s, params)`), and `r.bismark_args?.dedup` works when the caller omits it.
- **Input names can be reassigned in the script** (`args = args ?: ''`), so modules need no renamed locals. A local with the same name as an input (`def strandedness = ...`) is a lint error.
- **nf-schema works unchanged.** New optional samplesheet columns arrive as extra positional values from `samplesheetToList` (empty → `[]`), and an `object` param with `additionalProperties: {type: string}` validates `--args.<tool>=...`.
- **Remaining `meta`** is only nf-schema's transient meta map in the samplesheet parser, immediately flattened.
- **What users lose**: `ext.args` for tools the pipeline never configured (e.g. samtools stats/index) is no longer settable, and `withName` glob targeting is replaced by fixed tool keys. Per-sample columns need a schema entry per tool.
- **Verification gotcha**: `NXF_CACHE_DIR` is set globally here, so an isolated `NXF_HOME` does not isolate the cache; concurrent runs collide on the session lock. Run sequentially, or set `NXF_CACHE_DIR`/`NXF_WORK` per run.

## Not covered

- nf-test module/subworkflow tests and snapshots.
- GPU (parabricks) paths only linted, never run.
- Per-sample columns for tools other than trimgalore/bismark align (need schema entries).
