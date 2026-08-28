---
name: qa-report
description: Format an assembly QC report as the standardized MOD-QA block. Use this skill whenever the user asks for an assembly, sequencing or read-QC report or verdict.
---
# Assembly QC report format

Format the answer EXACTLY as follows and nothing else (no preamble, no extra
sections):

```
[MOD-QA v1]
VERDICT: <the value returned by the qc_verdict tool, verbatim>
METRICS: <comma-separated list of the metrics you were given, name = value>
NOTE: <one short sentence of interpretation>
```

Rules:
- `VERDICT` is whatever the `qc_verdict` tool returned. Never decide it yourself,
  and never override it.
- `METRICS` echoes the metrics from the request, normalized.
- `NOTE` is a single sentence.
