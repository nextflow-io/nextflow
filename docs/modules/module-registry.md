(module-registry-page)=

# Nextflow module registry

:::{versionadded} 26.04.0
:::

The [Nextflow module registry](https://registry.nextflow.io) is a central repository for sharing Nextflow modules.
It provides module discovery, version tracking, and integrity verification.

## Publishing modules

Use the `module publish` command to upload a module to the registry:

```console
$ nextflow module publish myorg/my-module
```

Provide a `namespace/name` reference (for an already-installed module) or a local directory path containing the module files.

The module directory must include the following files:

- `main.nf`: The module script.
- `meta.yml`: The module spec with name, version, description, and license.
- `README.md`: Module documentation.

Use `-dry-run` to validate the module structure without uploading:

```console
$ nextflow module publish myorg/my-module -dry-run
```

Use the {ref}`module spec <cli-module-spec>` and {ref}`module validate <cli-module-validate>` commands to prepare and verify a module before publishing.
See {ref}`dev-modules-page` for a guide on creating modules.

See {ref}`cli-module-publish` for the full command reference.

## Authentication

Publishing requires authentication.
Provide your registry token in one of two ways:

**Environment variable:**

```console
$ export NXF_REGISTRY_TOKEN=your-token
$ nextflow module publish myorg/my-module
```

**Nextflow configuration:**

```groovy
registry {
    apiKey = 'YOUR_REGISTRY_TOKEN'
}
```

## Registry configuration

By default, Nextflow uses the public registry at `https://registry.nextflow.io`.
Configure alternative or additional registries in your Nextflow configuration:

```groovy
registry {
    url = [
        'https://private.registry.myorg.com',
        'https://registry.nextflow.io'
    ]
    apiKey = '${MYORG_TOKEN}'
}
```

Nextflow queries registries in the order you specify until it finds a module.
The `apiKey` is used only for the primary (first) registry.

Specify a registry when publishing:

```console
$ nextflow module publish myorg/my-module -registry 'https://private.registry.myorg.com'
```
