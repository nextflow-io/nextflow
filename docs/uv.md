(uv-page)=

# uv environments

[uv](https://docs.astral.sh/uv/) is an extremely fast Python package and project manager, written in Rust. It can install Python packages, manage virtual environments, and handle Python versions.

Nextflow has built-in support for uv that allows the configuration of workflow dependencies using Python packages, requirements files, or pyproject.toml files.

This allows Nextflow applications to use Python packages managed by uv, taking advantage of its speed and reliability for creating reproducible Python environments.

## Prerequisites

This feature requires the [uv](https://docs.astral.sh/uv/getting-started/installation/) package manager to be installed on your system.

## How it works

Nextflow automatically creates and activates uv virtual environments given the dependencies specified by each process.

Dependencies are specified by using the {ref}`process-uv` directive, providing either the names of the required Python packages, the path of a requirements file, the path of a pyproject.toml file, or the path of an existing virtual environment directory.

You can specify the directory where the uv environments are stored using the `uv.cacheDir` configuration property (see the {ref}`configuration page <config-uv>` for details). When using a computing cluster, make sure to use a shared file system path accessible from all compute nodes.

:::{warning}
The uv environment feature is not supported by executors that use remote object storage as the work directory, e.g. AWS Batch.
:::

### Enabling uv environments

The use of uv packages specified using the {ref}`process-uv` directive needs to be enabled explicitly by setting the option shown below in the pipeline configuration file (i.e. `nextflow.config`):

```groovy
uv.enabled = true
```

Alternatively, it can be specified by setting the variable `NXF_UV_ENABLED=true` in your environment or by using the `-with-uv` command line option.

### Use Python package names

Python package names can be specified using the `uv` directive. Multiple package names can be specified by separating them with a blank space. For example:

```nextflow
process hello {
  uv 'numpy pandas matplotlib'

  script:
  '''
  python my_script.py
  '''
}
```

Using the above definition, a uv virtual environment that includes NumPy, Pandas, and Matplotlib is created and activated when the process is executed.

The usual pip package syntax and naming conventions can be used. The version of a package can be specified using pip version specifiers like so: `numpy>=1.24 pandas==2.0.0`.

### Use requirements files

uv environments can also be defined using a requirements file. For example, given a `requirements.txt` file:

```
numpy>=1.24.0
pandas>=2.0
scikit-learn
matplotlib
```

The environment for a process can be specified like so:

```nextflow
process hello {
  uv '/path/to/requirements.txt'

  script:
  '''
  python my_script.py
  '''
}
```

### Use pyproject.toml files

uv can also install dependencies from a `pyproject.toml` file:

```nextflow
process hello {
  uv '/path/to/pyproject.toml'

  script:
  '''
  python my_script.py
  '''
}
```

### Use existing environments

If you already have a uv virtual environment, you can use it directly by specifying the path:

```nextflow
process hello {
  uv '/path/to/existing/venv'

  script:
  '''
  python my_script.py
  '''
}
```

### Environment caching

Nextflow caches uv environments so that they are created only once for each unique set of packages. The cache directory can be configured using the `uv.cacheDir` setting or the `NXF_UV_CACHEDIR` environment variable.

### Python version

You can specify the Python version to use when creating virtual environments:

```groovy
uv.pythonVersion = '3.12'
```

### Advanced settings

The following settings are available in the `uv` scope of the Nextflow configuration:

- `uv.enabled`: Enable the use of uv environments (default: `false`)
- `uv.cacheDir`: The path where uv environments are stored
- `uv.createTimeout`: Timeout for environment creation (default: `20 min`)
- `uv.installOptions`: Extra command line options for `uv pip install`
- `uv.pythonVersion`: Python version for virtual environment creation
