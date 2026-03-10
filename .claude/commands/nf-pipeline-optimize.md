---
description: Analyze Nextflow pipelines for performance optimization opportunities, focusing on data flow patterns and pipeline architecture.
---

## User Input

```text
$ARGUMENTS
```

You **MUST** consider the user input before proceeding (if not empty). The user may provide:
- A path to a Nextflow pipeline repository or directory
- A specific `.nf` file to analyze
- Context about their execution environment (cloud, HPC, local)

## Goal

Analyze a Nextflow pipeline to identify performance optimization opportunities. This skill focuses on **data flow and architecture first** - where 90% of easy wins are found - rather than resource tuning (CPU/memory), which general-purpose AI tools already handle well.

**Key Insight**: Nextflow optimization is about **data flow first, resource tuning last**. Most pipeline slowdowns come from:
- Reading data multiple times
- Unnecessary serialization points
- Poor channel operator choices
- Suboptimal scatter/gather patterns

## Execution Steps

### 1. Locate Pipeline Files

First, identify the pipeline structure:

```bash
# Find all Nextflow files
find . -name "*.nf" -o -name "*.config" 2>/dev/null | head -50

# Check for standard pipeline structure
ls -la main.nf nextflow.config modules/ subworkflows/ 2>/dev/null || true
```

If the user provided a specific path, use that. Otherwise, search the current directory.

### 2. Read Pipeline Code

Read the main workflow files to understand the pipeline structure:
- `main.nf` or equivalent entry point
- Key module/process definitions
- Workflow definitions showing the DAG structure

### 3. Analyze Data Flow Anti-Patterns

**This is the highest-priority optimization category.** Look for:

#### 3.1. Redundant File Reads
Detect when the same input files are read by multiple processes independently:

```groovy
// ANTI-PATTERN: Reading the same BAM file multiple times
PROCESS_A(bam_ch)
PROCESS_B(bam_ch)
PROCESS_C(bam_ch)
// Each process independently reads the same file from storage

// BETTER: Read once, distribute via channel splitting
bam_ch
    .multiMap { bam ->
        for_a: bam
        for_b: bam
        for_c: bam
    }
    .set { split_bam }
PROCESS_A(split_bam.for_a)
PROCESS_B(split_bam.for_b)
PROCESS_C(split_bam.for_c)
```

**Check for**: Same channel variable used as input to multiple processes, especially for large files (BAM, FASTQ, VCF, CRAM).

#### 3.2. Missing Channel Splitting
Identify where `multiMap` or `branch` could distribute work more efficiently:

```groovy
// ANTI-PATTERN: Serial processing
results = ALIGN(reads)
CALL_VARIANTS(results)
CALCULATE_METRICS(results)
// No parallelism between independent operations

// BETTER: Parallel processing with multiMap
ALIGN(reads)
    .multiMap { aligned ->
        variants: aligned
        metrics: aligned
    }
    .set { aligned_split }
CALL_VARIANTS(aligned_split.variants)
CALCULATE_METRICS(aligned_split.metrics)
```

#### 3.3. Processes That Could Share Upstream Results
Flag when similar filtering/preparation is duplicated across processes:

```groovy
// ANTI-PATTERN: Each process does its own filtering
PROCESS_A(raw_data.filter { it.size() > 1000 })
PROCESS_B(raw_data.filter { it.size() > 1000 })

// BETTER: Filter once, share result
filtered = raw_data.filter { it.size() > 1000 }
PROCESS_A(filtered)
PROCESS_B(filtered)
```

### 4. Analyze Pipeline Architecture

#### 4.1. Bottleneck Detection (DAG Analysis)
Identify points where parallelism breaks down:

```groovy
// ANTI-PATTERN: Single-process bottleneck
parallel_results = MANY_PARALLEL_TASKS(inputs)
SINGLE_BOTTLENECK(parallel_results.collect())  // Everything waits here
MORE_PARALLEL(SINGLE_BOTTLENECK.out)
```

Look for processes that:
- Take `collect()` output and produce single output
- Have `maxForks = 1` unnecessarily
- Are the only path between parallel stages

#### 4.2. Premature Collection
Detect `collect()` or `groupTuple()` that serialize work too early:

```groovy
// ANTI-PATTERN: Premature collect
channel_of_files
    .collect()  // Waits for ALL files before proceeding
    | PROCESS_THAT_WORKS_ON_EACH

// BETTER: Process individually, collect only when needed
channel_of_files
    | PROCESS_THAT_WORKS_ON_EACH
    | collect  // Only collect for final aggregation
```

**Key patterns to flag**:
- `collect()` followed by `flatten()` (often unnecessary)
- `groupTuple()` with high cardinality keys early in pipeline
- Collection operators before processes that could work on individual items

#### 4.3. Missing Scatter-Gather Patterns
Identify processes that could be parallelized but aren't:

```groovy
// ANTI-PATTERN: Sequential processing of a list
process PROCESS_ALL {
    input: path(all_files)
    script: "for f in ${all_files}; do process $f; done"
}

// BETTER: Let Nextflow parallelize
process PROCESS_ONE {
    input: path(single_file)
    script: "process ${single_file}"
}
Channel.fromPath('data/*') | PROCESS_ONE
```

#### 4.4. Sequential Anti-Patterns
Flag chains that could use `branch` or parallel channels:

```groovy
// ANTI-PATTERN: Sequential decisions
if (condition) {
    PROCESS_A(data)
} else {
    PROCESS_B(data)
}

// BETTER: Use branch operator for parallel evaluation
data.branch {
    condition_a: it.property == 'a'
    condition_b: true  // fallback
}
.set { branched }
PROCESS_A(branched.condition_a)
PROCESS_B(branched.condition_b)
```

### 5. Analyze Scattering Efficiency

**Trade-off**: Too much splitting adds latency and overheads; too little causes long process times.

#### 5.1. Over-Scattering
Detect when tasks are too fine-grained:

```groovy
// ANTI-PATTERN: Scatter per-read (millions of tiny tasks)
reads.splitFastq(by: 1) | PROCESS_SINGLE_READ

// BETTER: Scatter in reasonable chunks
reads.splitFastq(by: 100000) | PROCESS_CHUNK
```

**Indicators of over-scattering**:
- `splitFastq`, `splitFasta`, `splitText` with very small `by:` values
- Process runtime dominated by staging/unstaging overhead
- Cloud execution with many tiny tasks (container startup overhead)

#### 5.2. Under-Scattering
Detect when parallelism opportunities are missed:

```groovy
// ANTI-PATTERN: Single-threaded chromosome processing
PROCESS_ALL_CHROMOSOMES(reference)

// BETTER: Scatter by chromosome
Channel.of(1..22, 'X', 'Y') | PROCESS_PER_CHROMOSOME
```

**Indicators of under-scattering**:
- Long-running processes that work on multiple independent units
- Processes that internally loop over items that could be parallelized

### 6. Analyze Environment-Specific Patterns

**Context matters**: Optimal patterns differ between cloud, HPC, and local execution.

#### 6.1. Compression/Decompression Overhead
Check for unnecessary compress/uncompress cycles:

```groovy
// ANTI-PATTERN on cloud (storage is fast, CPU is the bottleneck):
COMPRESS(output)   // gzip every intermediate
DECOMPRESS | NEXT_STEP

// BETTER on cloud: Skip compression for intermediates
// Only compress for long-term storage or final outputs
```

**Environment considerations**:
- **Cloud**: Storage I/O is often faster than compression CPU cost
- **HPC**: Compression may help with quota limits and slow filesystems
- **Local**: Usually skip compression for speed

#### 6.2. Combined Processes (Bash Piping)
Identify opportunities to combine I/O-heavy steps:

```groovy
// ANTI-PATTERN: Separate processes with intermediate files
SORT(bam)
INDEX(SORT.out)
FILTER(INDEX.out)

// BETTER: Combined with piping (harder to debug but faster)
process SORT_INDEX_FILTER {
    script:
    """
    samtools sort ${bam} | \
    samtools view -F 4 - | \
    samtools index - ${prefix}.sorted.filtered.bam
    """
}
```

**Trade-offs to note**:
- Combined processes are harder to debug and less reusable
- But they eliminate intermediate file I/O
- Best for stable, well-tested pipelines

### 7. Analyze Tool and Algorithm Choices

#### 7.1. Older Tools with Newer Alternatives
Flag commonly-known tool upgrades:

| Older Tool | Newer Alternative | Speedup Notes |
|------------|------------------|---------------|
| BWA MEM | BWA-MEM2, minimap2 | 2-3x faster alignment |
| STAR (standard) | STAR with shared genome | Memory efficiency |
| bowtie2 | bowtie2 --mm | Memory-mapped for multi-sample |
| gzip | pigz | Parallel compression |
| sort (coreutils) | GNU parallel sort | Multi-threaded |

#### 7.2. Algorithm Upgrades
Suggest newer algorithms where applicable:

| Traditional | Modern Alternative | Use Case |
|-------------|-------------------|----------|
| Full alignment | Pseudoalignment (kallisto, salmon) | RNA-seq quantification |
| Per-read alignment | Minimizer-based (minimap2) | Long reads |
| Exhaustive variant calling | ML-based (DeepVariant) | Germline calling |

**Note**: These are trade-offs - newer isn't always better for every use case.

### 8. Generate Optimization Report

Produce a structured report with findings organized by impact:

```markdown
# Pipeline Optimization Report

## Executive Summary
- Pipeline: [name]
- Total processes: [N]
- Critical issues: [N]
- Estimated improvement potential: [High/Medium/Low]

## Critical Issues (Fix These First)

### 1. Data Flow Anti-Patterns
| Issue | Location | Impact | Recommendation |
|-------|----------|--------|----------------|
| BAM read 5 times | main.nf:L42-78 | High I/O | Use multiMap to split channel |

### 2. Architecture Bottlenecks
| Issue | Location | Impact | Recommendation |
|-------|----------|--------|----------------|

## Medium Priority

### 3. Scattering Efficiency
| Issue | Location | Current | Suggested |
|-------|----------|---------|-----------|

### 4. Environment-Specific Optimizations
(Based on: [cloud/HPC/local])

## Low Priority / Consider

### 5. Tool and Algorithm Updates
| Current | Alternative | Trade-off |
|---------|-------------|-----------|

## Next Steps

1. [Most impactful change]
2. [Second most impactful]
3. [Third most impactful]
```

### 9. Provide Actionable Recommendations

For each finding:
1. Quote the specific code location
2. Explain WHY it's suboptimal
3. Show a concrete code fix
4. Note any trade-offs or considerations

## Important Notes

### What This Skill Does NOT Check
The following were intentionally excluded as they provide low signal-to-noise:

- **Generic resource tuning** (CPU/memory) - General AI handles this well
- **`storeDir` suggestions** - Often causes more problems than it solves
- **`scratch` directive** - Highly environment-specific, easy to misconfigure
- **`stageInMode` optimizations** - Rarely the bottleneck
- **Cache/resume debugging** - Better handled by Nextflow's built-in tools
- **Container pull optimization** - Minor impact compared to data flow issues

### Execution Environment Context
If the user specifies their environment:
- **Cloud (AWS/GCP/Azure)**: Prioritize reducing task count, skip compression
- **HPC (Slurm/PBS)**: Consider filesystem pressure, scheduler overhead
- **Local**: Focus on memory efficiency and disk I/O

### Trade-off Awareness
Always acknowledge trade-offs in recommendations:
- Combined processes save I/O but reduce debuggability
- Pseudoalignment is faster but loses positional information
- More scattering increases parallelism but adds overhead

## Context

$ARGUMENTS
