# Dataflow benchmarks

Whole-run benchmarks for the dataflow engine: wall time, peak RSS, and live heap. Rerun them after each phase of the GPars replacement and compare against the baseline.

| Script | Experiment | Workload |
|---|---|---|
| `dataflow.nf` | E1 | `n` items through `map`, `filter`, `flatMap`, `map`, `reduce` (default 10⁶) |
| `processes.nf` | E1 | one native (`exec:`) task per item (default 10⁴) |
| `slow-consumer.nf` | E2 | `n` three-field records into a `maxForks 1` process with a tiny native task (default 10⁶) |
| `slow-prepare.nf` | E3 | one native task per item whose `when` guard blocks for `ms` milliseconds (default 10³ tasks, 20 ms), next to a pure dataflow pipeline of `side` items (`side_ms`) |

## Running

```bash
make compile
modules/nf-dataflow/src/bench/run.sh -o /path/to/results                          # default matrix
modules/nf-dataflow/src/bench/run.sh -o /path/to/results dataflow:100000:256m:time  # one case
```

A case is `script:n:heap:passes`:
- `heap` is `default` or an `-Xmx` value.
- The `time` pass runs without sampling.
- The `heap` pass also runs `jcmd GC.class_histogram` every `SAMPLE` seconds (default 2) and keeps the largest total. Each sample forces a full GC, so the wall time of a `heap` pass is perturbed. Use the `time` pass for wall time.

Every run:
- goes through `launch.sh`, so it uses C1 and SerialGC like a development run;
- gets a fresh `NXF_HOME`, `NXF_CACHE_DIR` and work dir under `mktemp -d`, so it neither resumes nor contends for the user's cache;
- adds `-XX:+ExitOnOutOfMemoryError` and `-Xlog:gc` to `NXF_OPTS`, so an OOM ends the run instead of hanging it;
- is killed after `TIMEOUT` seconds (default 600).

Runs are sequential. Do not start another Nextflow run or a Gradle build while the benchmark runs.

## Output

`results.tsv` holds one row per pass:
- exit status; `oom` (`yes`, `no` or `timeout`);
- wall time and peak RSS (`/usr/bin/time -v`);
- peak live heap (largest full-GC histogram total);
- GC count, full GC count, total GC pause, and max heap after GC (from the GC log);
- the result line printed by the script.

Each pass directory holds the exact command, stdout, the GC log, the heap samples, and the class histogram at the peak. `machine.txt` records the host, CPU, memory, commit and `java -version`.
