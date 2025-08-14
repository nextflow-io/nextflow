(package-page)=

# Package Management

:::{versionadded} 25.04.0-edge
:::

Nextflow provides a unified package management system that allows you to specify dependencies using different package managers through a single, consistent interface. This system supports conda, pixi, and other package managers through a plugin-based architecture.

## Prerequisites

The unified package management system requires:
- The `preview.package` feature flag to be enabled
- The appropriate package manager installed on your system (conda, pixi, etc.)
- The corresponding Nextflow plugin for your chosen package manager

## How it works

Nextflow creates and activates the appropriate environment based on the package specifications and provider you choose. The system abstracts away the differences between package managers, providing a consistent interface regardless of the underlying tool.

## Enabling Package Management

The unified package management system is enabled using the `preview.package` feature flag:

```groovy
// nextflow.config
nextflow.preview.package = true
```

Alternatively, you can enable it with an environment variable:

```bash
export NXF_PREVIEW_PACKAGE=true
```

Or using a command-line option when running Nextflow:

```bash
nextflow run workflow.nf -c <(echo 'nextflow.preview.package = true')
```

## Basic Usage

### Package Directive

Use the `package` directive in your process definitions to specify dependencies:

```nextflow
process example {
    package "samtools=1.15 bcftools=1.15", provider: "conda"
    
    script:
    """
    samtools --help
    bcftools --help
    """
}
```

### Syntax

The basic syntax for the package directive is:

```nextflow
package "<package_specification>", provider: "<provider_name>"
```

- `<package_specification>`: Space-separated list of packages with optional version constraints
- `<provider_name>`: The package manager to use (e.g., "conda", "pixi")

### Multiple Packages

You can specify multiple packages in a single directive:

```nextflow
process analysis {
    package "bwa=0.7.17 samtools=1.15 bcftools=1.15", provider: "conda"
    
    script:
    """
    bwa mem ref.fa reads.fq | samtools view -bS - | bcftools view
    """
}
```

### Version Constraints

Different package managers support different version constraint syntaxes:

**Conda:**
```nextflow
package "python=3.9 numpy>=1.20 pandas<2.0", provider: "conda"
```

**Pixi:**
```nextflow
package "python=3.9 numpy>=1.20 pandas<2.0", provider: "pixi"
```

## Configuration

### Default Provider

You can set a default provider in your configuration:

```groovy
// nextflow.config
packages {
    provider = "conda"  // Default provider for all package directives
}
```

### Provider-Specific Settings

Each provider can have its own configuration:

```groovy
// nextflow.config
conda {
    enabled = true
    cacheDir = "$HOME/.nextflow/conda"
    channels = ['conda-forge', 'bioconda']
}

packages {
    provider = "conda"
    conda {
        channels = ['conda-forge', 'bioconda', 'defaults']
        useMicromamba = true
    }
    pixi {
        channels = ['conda-forge', 'bioconda']
    }
}
```

## Advanced Usage

### Environment Files

You can specify environment files instead of package lists:

```nextflow
process fromFile {
    package file("environment.yml"), provider: "conda"
    
    script:
    """
    python analysis.py
    """
}
```

### Per-Provider Options

Some providers support additional options:

```nextflow
process withOptions {
    package "biopython scikit-learn", 
           provider: "conda", 
           channels: ["conda-forge", "bioconda"]
    
    script:
    """
    python -c "import Bio; import sklearn"
    """
}
```

## Supported Providers

### Conda

The conda provider supports:
- Package specifications with version constraints
- Custom channels
- Environment files (`.yml`, `.yaml`)
- Micromamba as an alternative backend

```nextflow
process condaExample {
    package "bioconda::samtools=1.15 conda-forge::numpy", 
           provider: "conda"
    
    script:
    """
    samtools --version
    python -c "import numpy; print(numpy.__version__)"
    """
}
```

### Pixi

The pixi provider supports:
- Package specifications compatible with pixi
- Custom channels
- Project-based environments

```nextflow
process pixiExample {
    package "samtools bcftools", provider: "pixi"
    
    script:
    """
    samtools --version
    bcftools --version
    """
}
```

## Migration from Legacy Directives

### From conda directive

**Before (legacy):**
```nextflow
process oldWay {
    conda "samtools=1.15 bcftools=1.15"
    
    script:
    "samtools --help"
}
```

**After (unified):**
```nextflow
process newWay {
    package "samtools=1.15 bcftools=1.15", provider: "conda"
    
    script:
    "samtools --help"
}
```

### Deprecation Warnings

When the unified package management system is enabled, using the legacy `conda` directive will show a deprecation warning:

```
WARN: The 'conda' directive is deprecated when preview.package is enabled. 
      Use 'package "samtools=1.15", provider: "conda"' instead
```

## Best Practices

### 1. Pin Package Versions

Always specify exact versions for reproducibility:

```nextflow
// Good
package "samtools=1.15 bcftools=1.15", provider: "conda"

// Avoid (may cause reproducibility issues)
package "samtools bcftools", provider: "conda"
```

### 2. Use Appropriate Channels

Specify the most appropriate channels for your packages:

```nextflow
process bioinformatics {
    package "bioconda::samtools conda-forge::pandas", provider: "conda"
    
    script:
    """
    samtools --version
    python -c "import pandas"
    """
}
```

### 3. Group Related Packages

Keep related packages together in the same environment:

```nextflow
process genomicsAnalysis {
    package "samtools=1.15 bcftools=1.15 htslib=1.15", provider: "conda"
    
    script:
    """
    # All tools are from the same suite and work well together
    samtools view input.bam | bcftools view
    """
}
```

### 4. Test Your Environments

Always test your package environments before deploying:

```bash
# Test package resolution
nextflow run test.nf --dry-run -preview

# Test actual execution
nextflow run test.nf -resume
```

## Troubleshooting

### Common Issues

**Package not found:**
- Check package name spelling
- Verify the package exists in specified channels
- Try different channels or provider

**Version conflicts:**
- Relax version constraints if possible
- Check for incompatible package combinations
- Consider using a different provider

**Slow environment creation:**
- Use `useMicromamba = true` for faster conda operations
- Consider pre-built environments
- Use appropriate cache directories

### Environment Inspection

You can inspect created environments using provider-specific commands:

```bash
# For conda environments
conda env list
conda list -n nextflow-env-hash

# For pixi environments
pixi info
```

## Integration with Wave

The package management system integrates seamlessly with Wave containers. When Wave is enabled, environments are automatically containerized:

```groovy
// nextflow.config
wave.enabled = true
nextflow.preview.package = true
```

```nextflow
process containerized {
    package "samtools=1.15", provider: "conda"
    
    script:
    """
    # This runs in a Wave container with samtools pre-installed
    samtools --version
    """
}
```

## Limitations

- The unified package management system is currently in preview
- Plugin availability may vary for different providers
- Some legacy features may not be fully supported yet
- Provider-specific options may be limited

## See Also

- {ref}`conda-page` - Legacy conda directive documentation
- {ref}`config-packages` - Package management configuration options
- {ref}`wave-page` - Wave container integration