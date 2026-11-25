#!/usr/bin/env python3

"""Contract fixture for native external Python scripts.

This file should lint/type-check as ordinary Python. The only Nextflow-specific
surface is the importable runtime object, backed by the task sidecar generated
by the engine.
"""

from __future__ import annotations

import json
from pathlib import Path

from nextflow.script import nextflow


def count_records(fasta: Path) -> int:
    return sum(1 for line in fasta.read_text().splitlines() if line.startswith(">"))


reads = Path(nextflow.input["reads"])
summary = Path(nextflow.output["summary"])
context = Path(nextflow.output["context"])
prefix = nextflow.params["prefix"]

summary.write_text(f"{prefix}\trecords\t{count_records(reads)}\n")
context.write_text(
    json.dumps(
        {
            "language": "python",
            "reads": str(reads),
            "summary": str(summary),
            "prefix": prefix,
            "cpus": nextflow.cpus,
            "attempt": nextflow.attempt,
        },
        indent=2,
        sort_keys=True,
    )
    + "\n"
)
