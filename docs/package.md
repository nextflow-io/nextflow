(package-page)=

# Package Management

:::{versionadded} 25.07.0-edge
:::

Nextflow provides a unified package management system that allows you to specify dependencies using different package managers through a single, consistent interface. This system supports conda, pixi, uv, and other package managers through a plugin-based architecture.

## Prerequisites

The unified package management system requires:
- The `preview.package` feature flag to be enabled
- The appropriate package manager installed on your system (conda, pixi, uv, etc.)
- The corresponding Nextflow plugin for your chosen package manager (`nf-conda`, `nf-pixi`, `nf-uv`, etc.)

## How it works

Nextflow creates and activates the appropriate environment based on the package specifications and provider you choose. The system abstracts away the differences between package managers, providing a consistent interface regardless of the underlying tool.

## Enabling Package Management

:::{important}
The `package` directive requires the **v2 syntax parser**, because `package` is a reserved word in the legacy (v1) parser. Enable it by setting `export NXF_SYNTAX_PARSER=v2` before running Nextflow.
:::

The unified package management system is enabled using the `nextflow.preview.package` feature flag, either at the top of your pipeline script:

```nextflow
// main.nf
nextflow.preview.package = true
```

Or in your configuration file:

```groovy
// nextflow.config
nextflow.preview.package = true
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

### Automatic environment file detection

When a process does not declare a `package` directive, Nextflow automatically looks in the process **module directory** for a known manifest file and uses it, analogous to how Wave automatically picks up a `Dockerfile`. This means you can drop a manifest file next to your module and omit the directive entirely:

```
modules/
└── myprocess/
    ├── main.nf            # process with no `package` directive
    └── environment.yml    # auto-detected -> provider: conda
```

The first manifest found (providers are scanned in name order) determines the provider:

| Manifest file | Provider | Resolved with |
|---|---|---|
| `environment.yml`, `environment.yaml` | `conda` | `conda env create -f` |
| `pixi.toml` | `pixi` | pixi project install |
| `requirements.txt`, `pyproject.toml` | `uv` | `uv pip install -r` |
| `manifest.scm`, `guix.scm` | `guix` | `guix package --manifest=` |
| `DESCRIPTION` | `pak` | `pak::local_install_deps()` |

`nix` and `install2r` do not participate in auto-detection: a `flake.nix`/`shell.nix` can describe a dev shell or arbitrary packages (no unambiguous "install these tools" mapping), and `install2.r` is a CLI for named packages with no manifest-file convention. Use an explicit `package` directive for those.

Auto-detection is enabled by default when the package preview feature is on. Disable it with:

```groovy
// nextflow.config
packages {
    autoDetect = false
}
```

An explicit `package` directive always takes precedence over auto-detection. Other providers expose their manifest file names through the provider plugin API (`getManifestFileNames()`), so additional managers can opt in to auto-detection.

### Custom-named manifest files

Auto-detection only looks for the conventional file names above. If your manifest has a different name (e.g. `env.yml` instead of `environment.yml`), pass an explicit path to the directive and the provider will still recognise it, subject to the file-type rules each tool understands:

```nextflow
process customManifest {
    package "${moduleDir}/env.yml", provider: "conda"   // custom name, still a conda YAML

    script:
    """
    python analysis.py
    """
}
```

| Provider | Recognised as a manifest file when the path ends with |
|---|---|
| `conda` | `.yml`, `.yaml` (conda env file) or `.txt` (package list), any base name; a remote `.yml`/`.yaml` URL also works |
| `uv` | `.txt`, `.in` (requirements), any base name; or a path ending in `pyproject.toml` |
| `pixi` | `.toml`, `.lock`, any base name |
| `guix` | any existing file (used as a `--manifest` file) |
| `pak` | any existing file (installed as a package `DESCRIPTION` in its directory) |
| `nix`, `install2r` | manifest files are **not** supported; pass explicit package/flake refs or names |

A path that doesn't match these rules is treated as a literal package specification.

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

### Per-process options (override config)

Each provider reads its install/create options from its config scope by default
(e.g. `uv.installOptions`, `conda.createOptions`). A process can override that
default for its own environment using the `options` map in the directive. The
key matches the provider's config option name (`installOptions` for uv, nix,
guix, pak and install2r; `createOptions` for conda and pixi), so it's "the same
knob, two scopes":

```nextflow
process fastResolve {
    // overrides uv.installOptions for this process only; config is the fallback
    package "numpy pandas", provider: "uv", options: [installOptions: "--no-cache"]

    script:
    """
    python -c "import numpy, pandas"
    """
}
```

A per-process override participates in the environment cache key, so two
processes that request the same packages with different options get distinct
environments.

## Supported Providers

Each provider is supplied by its own `nf-<provider>` plugin. Load the plugins you need in your configuration (e.g. `plugins { id 'nf-uv' }`). The table below lists the supported package managers and whether they can be built into containers by [Wave](https://docs.seqera.io/wave); providers that are not Wave-compatible create environments on the local file system and are skipped by Wave so the provider plugin resolves them locally.

| Provider (`provider:`) | Plugin | Manages | Wave-compatible |
|---|---|---|---|
| `conda` | `nf-conda` | Conda packages / environment files | Yes |
| `mamba` / `micromamba` | `nf-conda` | Conda packages (mamba backend) | Yes |
| `pixi` | `nf-pixi` | Conda packages (pixi) | Yes (built as conda) |
| `uv` | `nf-uv` | Python packages (uv virtual environments) | No (local only) |
| `nix` | `nf-nix` | Nix packages / flake refs | No (local only) |
| `guix` | `nf-guix` | GNU Guix packages | No (local only) |
| `pak` | `nf-pak` | R packages (pak) | No (local only) |
| `install2r` | `nf-install2r` | R packages (install2.r / littler) | No (local only) |

:::{note}
When `wave.enabled` is set, only Wave-compatible providers are built into containers. Processes that use a local-only provider (uv, nix, guix, pak, install2r) are not sent to Wave; their environments are created on the local file system by the corresponding plugin instead, so conda-on-Wave and uv-locally can coexist in the same pipeline.
:::


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

### uv

The [uv](https://docs.astral.sh/uv/) provider manages Python virtual environments. It is provided by the `nf-uv` plugin and supports:
- Python package lists with pip-style version specifiers (e.g. `numpy>=1.24 pandas==2.0.0`)
- `requirements.txt` / `requirements.in` files
- `pyproject.toml` files
- Existing virtual environment directories

```nextflow
process uvExample {
    package "numpy pandas matplotlib", provider: "uv"

    script:
    """
    python -c "import numpy, pandas, matplotlib; print(numpy.__version__)"
    """
}
```

You can also point the directive at a requirements or project file:

```nextflow
process uvFromFile {
    package "/path/to/requirements.txt", provider: "uv"

    script:
    """
    python analysis.py
    """
}
```

The uv provider can be configured through the `uv` config scope:

```groovy
// nextflow.config
uv {
    cacheDir = "$HOME/.nextflow/uv"
    pythonVersion = '3.12'
    installOptions = '--no-cache'
    createTimeout = '20 min'
}
```

:::{note}
The uv provider creates environments on the local file system and is not supported by executors that use remote object storage as the work directory (e.g. AWS Batch). Use a POSIX-compatible work directory, or set `uv.cacheDir` to a shared file-system path accessible from all compute nodes.

The uv provider is also not supported by Wave container builds, since Wave only builds conda-based package environments. Use the uv provider with local or HPC executors rather than enabling `wave.enabled` together with `provider: "uv"`.
:::

### Nix

The [Nix](https://nixos.org/) provider manages dependencies as Nix profiles. It is provided by the `nf-nix` plugin and supports:
- Bare package names, resolved against the `nixpkgs` flake by default (e.g. `samtools` becomes `nixpkgs#samtools`)
- Fully-qualified flake references (e.g. `nixpkgs#hello`, `github:owner/repo#pkg`)

```nextflow
process nixExample {
    package "hello", provider: "nix"

    script:
    """
    hello
    """
}
```

The flake used to resolve bare names, and other settings, are configured through the `nix` config scope:

```groovy
// nextflow.config
nix {
    cacheDir = "$HOME/.nextflow/nix"
    flakeRef = 'nixpkgs'        // flake used to resolve bare package names
    installOptions = ''         // extra args for `nix profile install`
    createTimeout = '20 min'
}
```

The Nix provider requires the `nix` command with the `nix-command` and `flakes` experimental features (Nextflow passes these on the command line). Environments are created on the local file system and are not supported by Wave or remote object-storage work directories.

### Guix

The [GNU Guix](https://guix.gnu.org/) provider manages dependencies as Guix profiles. It is provided by the `nf-guix` plugin and supports package specifications resolved by `guix package`:

```nextflow
process guixExample {
    package "hello", provider: "guix"

    script:
    """
    hello
    """
}
```

The Guix provider is configured through the `guix` config scope:

```groovy
// nextflow.config
guix {
    cacheDir = "$HOME/.nextflow/guix"
    installOptions = ''         // extra args for `guix package --install`
    createTimeout = '20 min'
}
```

The Guix provider requires the `guix` command. Environments are created on the local file system and are not supported by Wave or remote object-storage work directories.

### R: pak

The [pak](https://pak.r-lib.org/) provider installs R packages into a per-environment R library using `pak::pkg_install()`. It is provided by the `nf-pak` plugin:

```nextflow
process pakExample {
    package "dplyr ggplot2", provider: "pak"

    script:
    """
    Rscript -e 'library(dplyr); library(ggplot2)'
    """
}
```

The pak provider activates the environment by setting `R_LIBS_USER` to the created library. It is configured through the `pak` config scope:

```groovy
// nextflow.config
pak {
    cacheDir = "$HOME/.nextflow/pak"
    installOptions = ''         // extra args appended to pak::pkg_install()
    createTimeout = '20 min'
}
```

The pak provider requires `R` (the `Rscript` command) with the `pak` package available. Environments are created on the local file system and are not supported by Wave or remote object-storage work directories.

### R: install2.r

The `install2.r` provider (from the [littler](https://github.com/eddelbuettel/littler) package) installs CRAN packages into a per-environment R library. It is provided by the `nf-install2r` plugin:

```nextflow
process install2rExample {
    package "dplyr ggplot2", provider: "install2r"

    script:
    """
    Rscript -e 'library(dplyr); library(ggplot2)'
    """
}
```

Like pak, it activates via `R_LIBS_USER` and is configured through the `install2r` config scope:

```groovy
// nextflow.config
install2r {
    cacheDir = "$HOME/.nextflow/install2r"
    installOptions = ''         // extra install2.r CLI flags, e.g. '--error'
    createTimeout = '20 min'
}
```

The install2r provider requires the `install2.r` command on the `PATH`. Environments are created on the local file system and are not supported by Wave or remote object-storage work directories.

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
# Preview the workflow without executing any task
nextflow run test.nf -preview

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

The package management system integrates with Wave containers. Wave can only build conda-based package environments, so when Wave is enabled, processes using the `conda`, `mamba`, `micromamba`, or `pixi` providers are automatically containerized, while processes using local-only providers (`uv`, `nix`, `guix`, `pak`, `install2r`) skip Wave and their environments are created locally by the corresponding provider plugin. This allows mixed pipelines, e.g. conda processes built by Wave alongside uv processes created locally:

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

## Plugin Architecture

The system is extensible through plugins. Each package manager is implemented as a plugin that provides a `PackageProviderExtension`:

- `nf-conda` - Conda support
- `nf-pixi` - Pixi support
- `nf-uv` - uv (Python) support
- `nf-nix` - Nix support
- `nf-guix` - GNU Guix support
- `nf-pak` - R pak support
- `nf-install2r` - R install2.r (littler) support

Custom package managers can be added by implementing the `PackageProvider` interface and registering it as a plugin.

## Limitations

- The unified package management system is currently in preview
- Plugin availability may vary for different providers
- Some legacy features may not be fully supported yet
- Provider-specific options may be limited

## See Also

- {ref}`conda-page` - Legacy conda directive documentation
- {ref}`config-packages` - Package management configuration options
- {ref}`wave-page` - Wave container integration