---
name: sequence-report
description: Format a sequencing or genome-assembly summary as a standardized QC report. Use this skill whenever the user asks for a sequence, assembly, or read-QC summary or verdict.
---
# Sequence QC report format

When producing a sequencing or assembly summary, format the answer EXACTLY as
follows and nothing else (no preamble, no extra sections):

```
[SEQ-REPORT v1]
STATUS: <PASS | WARN | FAIL>
METRICS: <comma-separated list of the key metrics you were given, e.g. N50, total length, GC%>
NOTE: <one short sentence of interpretation>
```

Rules:
- `STATUS` is your overall verdict from the metrics provided.
- `METRICS` echoes the metrics from the request, normalized (name = value).
- `NOTE` is a single terse sentence — no more.
