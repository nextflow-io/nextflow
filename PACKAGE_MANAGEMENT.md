# Unified Package Management System

This document describes the new unified package management system in Nextflow, introduced as a preview feature behind the `nextflow.preview.package` flag.

## Overview

The unified package management system provides a consistent interface for managing packages across different package managers (conda, pixi, mamba, etc.) through a single `package` directive.

## Enabling the Feature

Add this to your `nextflow.config`:

```groovy
nextflow.preview.package = true
```

## Basic Usage

### Single Package

```groovy
process example {
    package "samtools=1.17", provider: "conda"
    
    script:
    """
    samtools --version
    """
}
```

### Multiple Packages

```groovy
process example {
    package ["samtools=1.17", "bcftools=1.18"], provider: "conda"
    
    script:
    """
    samtools --version
    bcftools --version
    """
}
```

### Using Default Provider

Configure a default provider in your config:

```groovy
packages {
    provider = 'conda'
}
```

Then use:

```groovy
process example {
    package "samtools=1.17"  // uses default provider
    
    script:
    """
    samtools --version
    """
}
```

### Advanced Configuration

```groovy
process example {
    package {
        provider = "conda"
        packages = ["samtools=1.17", "bcftools=1.18"]
        channels = ["conda-forge", "bioconda"]
        options = [
            createTimeout: "30 min"
        ]
    }
    
    script:
    """
    samtools --version
    bcftools --version
    """
}
```

### Environment Files

```groovy
process example {
    package {
        provider = "conda"
        environment = file("environment.yml")
    }
    
    script:
    """
    python script.py
    """
}
```

## Supported Providers

- `conda` - Anaconda/Miniconda package manager
- `pixi` - Fast conda alternative with lockfiles
- `mamba` - Fast conda alternative
- `micromamba` - Minimal conda implementation

## Configuration

### Global Configuration

```groovy
// nextflow.config
nextflow.preview.package = true

packages {
    provider = 'conda'  // default provider
}

// Provider-specific configurations
conda {
    channels = ['conda-forge', 'bioconda']
    createTimeout = '20 min'
}

pixi {
    cacheDir = '/tmp/pixi-cache'
}
```

## Wave Integration

The unified package system integrates with Wave for containerization:

```groovy
process example {
    package "samtools=1.17", provider: "conda"
    
    script:
    """
    samtools --version
    """
}
```

Wave will automatically create a container with the specified packages.

## Backward Compatibility

Old `conda` and `pixi` directives continue to work but show deprecation warnings when the preview feature is enabled:

```groovy
process oldStyle {
    conda 'samtools=1.17'  // Shows deprecation warning
    
    script:
    """
    samtools --version
    """
}
```

## Migration Guide

### From conda directive

**Before:**
```groovy
process example {
    conda 'samtools=1.17 bcftools=1.18'
    script: "samtools --version"
}
```

**After:**
```groovy
process example {
    package ["samtools=1.17", "bcftools=1.18"], provider: "conda"
    script: "samtools --version"
}
```

### From pixi directive

**Before:**
```groovy
process example {
    pixi 'samtools bcftools'
    script: "samtools --version"
}
```

**After:**
```groovy
process example {
    package ["samtools", "bcftools"], provider: "pixi"
    script: "samtools --version"
}
```

## Plugin Architecture

The system is extensible through plugins. Package managers are implemented as plugins that extend the `PackageProviderExtension` interface:

- `nf-conda` - Conda support
- `nf-pixi` - Pixi support

Custom package managers can be added by implementing the `PackageProvider` interface and registering as a plugin.

## Examples

See the test files for complete examples:
- `tests/package-test.nf` - Basic usage examples
- `tests/integration-test.nf` - Integration and backward compatibility tests