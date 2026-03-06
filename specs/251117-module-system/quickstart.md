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

This downloads the module to `modules/@nf-core/fastqc/` and updates `nextflow_spec.json` with the installed version.

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

# Run specific version
nextflow module run nf-core/fastqc --input 'data/*.fastq.gz' -version 1.0.0

# With Nextflow options
nextflow module run nf-core/salmon \
    --reads reads.fq \
    --index salmon_index \
    -profile docker \
    -resume
```

## 3. View Module Information

```bash
# Show module metadata and a generated usage template
nextflow module info nf-core/fastqc

# Show a specific version
nextflow module info nf-core/fastqc -version 1.0.0

# JSON output for scripting
nextflow module info nf-core/fastqc -json
```

---

## 4. Manage Module Versions

### Version tracking

Module versions are automatically recorded in `nextflow_spec.json` by `nextflow module install`. You can also pin versions manually:

```json
// nextflow_spec.json
{
  "modules": {
    "@nf-core/fastqc": "1.0.0",
    "@nf-core/bwa-align": "1.2.0"
  }
}
```

Alternatively, declare versions in `nextflow.config` (not currently used):

```nextflow
modules {
    '@nf-core/fastqc' = '1.0.0'
    '@nf-core/bwa-align' = '1.2.0'
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

Change the version in `nextflow_spec.json` (or `nextflow.config`), then run your workflow. Nextflow automatically downloads the new version.

---

## 5. Search for Modules

```bash
# Search by keyword
nextflow module search alignment

# Limit results
nextflow module search "quality control" -limit 5

# JSON output for scripting
nextflow module search bwa -json
```

---

## 6. Work with Private Registries

### Configure authentication

```nextflow
// nextflow.config
registry {
    // Multiple registries (tried in order)
    url = [
        'https://private.registry.myorg.com',
        'https://registry.nextflow.io/api'
    ]
    apiKey = 'MYORG_TOKEN'  // Applied to the primary (first) registry only
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