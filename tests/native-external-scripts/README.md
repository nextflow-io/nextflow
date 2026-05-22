# Native external scripts fixture

This directory sketches the executable test shape for `adr/20260522-native-external-scripts.md`, following Snakemake's `tests/test_script` style.

It is intentionally **not wired into the current `tests/checks/run.sh` sweep** because the proposed `script: file '...'` runtime and context injection do not exist yet. Once implemented, the fixture can be promoted into the normal checks by invoking `main.nf` and comparing stdout with `checks/expected.txt`.

The important contract is the file layout and source shape:

- `main.nf` declares processes that point at external script assets.
- `scripts/bash_sort.sh` is valid Bash as written and expects `nextflow_*` associative arrays.
- `scripts/python_summary.py` is valid Python as written and imports `nextflow.script.nextflow`.
- `scripts/r_summary.R` is valid R as written and expects an injected `nextflow` object.
