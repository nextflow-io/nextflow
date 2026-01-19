# Quickstart: Nextflow Module System

This guide covers the essential workflows for using the Nextflow module system.

## Prerequisites

- Nextflow 25.x or later (with module system support)
- Network connectivity for initial module downloads
- Optional: `NXF_REGISTRY_TOKEN` for publishing

---

## 1. Install and Use a Module

### Install a module

```bash
# Install latest version
nextflow module install nf-core/fastqc

# Install specific version
nextflow module install nf-core/fastqc -version 1.0.0
```

This downloads the module to `modules/@nf-core/fastqc/` and updates `nextflow.config`.

### Use in your workflow

```groovy
// main.nf
include { FASTQC } from '@nf-core/fastqc'

workflow {
    reads = Channel.fromFilePairs('data/*_{1,2}.fastq.gz')
    FASTQC(reads)
}
```

### Run your workflow

```bash
nextflow run main.nf
```

---

## 2. Run a Module Directly

Execute a module without writing a wrapper workflow:

```bash
# Basic usage
nextflow module run nf-core/fastqc --input 'data/*.fastq.gz'

# With tool arguments
nextflow module run nf-core/bwa-align \
    --reads 'samples/*_{1,2}.fastq.gz' \
    --reference genome.fa \
    --tools:bwa:K 100000000

# With Nextflow options
nextflow module run nf-core/salmon \
    --reads reads.fq \
    --index salmon_index \
    -profile docker \
    -resume
```

---

## 3. Manage Module Versions

### Configure versions in nextflow.config

```groovy
// nextflow.config
modules {
    '@nf-core/fastqc' = '1.0.0'
    '@nf-core/bwa-align' = '1.2.0'
    '@nf-core/samtools' = '2.1.0'
}
```

### Check module status

```bash
# List all modules
nextflow module list

# Output:
# MODULE                  CONFIGURED  INSTALLED  LATEST  STATUS
# @nf-core/fastqc         1.0.0       1.0.0      1.2.0   outdated
# @nf-core/bwa-align      1.2.0       1.2.0      1.2.0   up-to-date
# @nf-core/samtools       2.1.0       -          2.1.0   missing
```

### Update a module

Change the version in `nextflow.config`, then run your workflow. Nextflow automatically downloads the new version.

```groovy
modules {
    '@nf-core/fastqc' = '1.2.0'  // Changed from 1.0.0
}
```

---

## 4. Search for Modules

```bash
# Search by keyword
nextflow module search alignment

# Limit results
nextflow module search "quality control" -limit 5

# JSON output for scripting
nextflow module search bwa -json
```

---

## 5. Configure Tool Arguments

### Define in meta.yaml (module author)

```yaml
# modules/@nf-core/bwa-align/meta.yaml
tools:
  - bwa:
      description: BWA aligner
      args:
        K:
          flag: "-K"
          type: integer
          description: "Process INT input bases in each batch"
        Y:
          flag: "-Y"
          type: boolean
          description: "Use soft clipping for supplementary alignments"
```

### Configure in nextflow.config (user)

```groovy
// nextflow.config
process {
    withName: 'BWA_ALIGN' {
        tools.bwa.args.K = 100000000
        tools.bwa.args.Y = true
    }
}
```

### Access in script (module author)

```groovy
// main.nf
process BWA_ALIGN {
    script:
    """
    bwa mem ${tools.bwa.args} -t $task.cpus $index $reads
    """
}
```

---

## 6. Work with Private Registries

### Configure authentication

```groovy
// nextflow.config
registry {
    // Multiple registries (tried in order)
    url = [
        'https://private.registry.myorg.com',
        'https://registry.nextflow.io'
    ]

    auth {
        'private.registry.myorg.com' = '${MYORG_TOKEN}'
        'registry.nextflow.io' = '${NXF_REGISTRY_TOKEN}'
    }
}
```

### Or use environment variable

```bash
export NXF_REGISTRY_TOKEN=your-token-here
nextflow module install nf-core/fastqc
```

---

## 7. Publish a Module

### Prepare your module

```
my-module/
├── main.nf          # Required: entry point
├── meta.yaml        # Required for registry
├── README.md        # Required for registry
└── tests/           # Recommended
```

### Validate before publishing

```bash
nextflow module publish myorg/my-module -dry-run
```

### Publish to registry

```bash
export NXF_REGISTRY_TOKEN=your-token
nextflow module publish myorg/my-module
```

---

## 8. Handle Local Modifications

If you modify a module locally (for debugging), Nextflow protects your changes:

```bash
# This warns and does NOT override your changes
nextflow module install nf-core/fastqc -version 1.1.0
# Warning: Module @nf-core/fastqc has local modifications. Use -force to override.

# Force replacement if needed
nextflow module install nf-core/fastqc -version 1.1.0 -force
```

---

## 9. Remove a Module

```bash
# Remove module and config entry
nextflow module remove nf-core/fastqc

# Keep config entry (just delete local files)
nextflow module remove nf-core/fastqc -keep-config

# Keep local files (just remove from config)
nextflow module remove nf-core/fastqc -keep-files
```

---

## Common Patterns

### Install all configured modules

```bash
# Installs all modules listed in nextflow.config
nextflow module install
```

### Offline operation

Modules are cached locally in `modules/`. Once installed, workflows run without network access.

### Git integration

The `modules/` directory is intended to be committed to your git repository:

```bash
git add modules/
git commit -m "Add module dependencies"
```

---

## Troubleshooting

### Module not found

```bash
# Check if module exists in registry
nextflow module search exact-module-name

# Verify spelling and scope
# Correct: @nf-core/fastqc
# Wrong:   @nfcore/fastqc, nf-core/fastqc (without @)
```

### Authentication errors

```bash
# Verify token is set
echo $NXF_REGISTRY_TOKEN

# Check registry config
grep -A5 'registry' nextflow.config
```

### Version conflicts

If two modules require incompatible versions of a dependency:
- Nextflow selects the highest compatible version automatically
- If no compatible version exists, an error lists the conflicts

### Checksum warnings

```
Warning: Module @nf-core/fastqc has local modifications
```

This means the local module content differs from the registry version. Your changes are preserved. Use `-force` only if you want to discard local changes.