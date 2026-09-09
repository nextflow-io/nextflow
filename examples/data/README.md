# Shared example fixtures

Input data shared by the examples, kept here once instead of per example. The examples that
need it reach it through a `data` symlink pointing at this directory, so all of them read the
same bytes and one fetch serves every example:

```
examples/agents/07_module-as-tool/data -> ../../data
examples/agents/09_goal-directed/data  -> ../../data
examples/agents/11_contig-filter/data  -> ../../data
examples/agents/12_isolate-triage/data -> ../../data
```

The files are not committed — see `.gitignore` for why. Fetch them once, from the repository
root:

```bash
mkdir -p examples/data
curl -sL https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/sarscov2/illumina/fastq/sarscov2_mus-musculus.fastq.gz \
  -o examples/data/sample.fastq.gz
gunzip -c examples/data/sample.fastq.gz > examples/data/sample.fastq
```

That produces both forms, and both are needed:

| File | Used by | Why this form |
|---|---|---|
| `sample.fastq` | `07_module-as-tool`, `09_goal-directed`, `12_isolate-triage` | uncompressed |
| `sample.fastq.gz` | `11_contig-filter` | `SEQTK_SAMPLE` names its output after its input and emits `*.fastq.gz`, so the input must be gzipped |

## What the data is

The sarscov2 Illumina read set from [nf-core/test-datasets](https://github.com/nf-core/test-datasets)
(`modules` branch), 136,166 single-end reads. It is a small **viral** read set, which is why
`09_goal-directed` and `12_isolate-triage` document a *failing* QC verdict as their expected
outcome: it assembles to roughly N50=310 over ~6 contigs, below the N50>=500 bar those examples
gate on. Substituting a different dataset invalidates the expected results recorded in those
examples and in `examples/agents/README.md`.

`19_shell-tools/` is the exception to all of this: it ships its own 2.9 kB FASTA in its own
`data/` directory, because the statistics quoted in its README only reproduce with those exact
bytes.
