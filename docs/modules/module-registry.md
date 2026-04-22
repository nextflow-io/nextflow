(module-registry-page)=

# Nextflow module registry

:::{versionadded} 26.04.0
:::

The [Nextflow module registry](https://registry.nextflow.io) is a central repository for sharing Nextflow modules.
It provides module discovery, version tracking, and integrity verification.

The module registry shares its login, access tokens, and web interface with the {ref}`Nextflow plugin registry <plugin-registry-page>`.

(module-registry-claim)=

## Claiming a namespace

A namespace is an organizational prefix (e.g., the namespace `nf-core` owns modules like `nf-core/fastqc` and `nf-core/samtools`).
As a namespace member, you can publish new modules and release new versions under that prefix.
Each namespace can have multiple members who share publishing rights.

To claim a namespace:

1. Open the [Nextflow module registry](https://registry.nextflow.io/) in a browser.

2. Log in to [Seqera](https://cloud.seqera.io/login) with your GitHub or Google account, or by providing an email address.

    :::{note}
    If you are logging in for the first time, Seqera sends an authentication link to your email address to complete the login process.
    :::

3. Go to **My Resources** and select the **My Namespaces** tab.

4. Select **Claim Namespace**

4. Enter your **Namespace name**, **Organization**, **Project URL**, and a **Description**.

5. Select **Submit Request**.

The namespace shows as **PENDING REVIEW** until an administrator approves the request.
Admin approval is required only once.
Once approved, you can publish modules using the `nextflow module publish` command.

(module-registry-access-token)=

## Creating an access token

An API access token is required to publish modules to the registry.

To create an API access token:

1. Open the [Nextflow module registry](https://registry.nextflow.io/) in a browser.

2. Log in to [Seqera](https://cloud.seqera.io/login) with your GitHub or Google account, or by providing an email address.

    :::{note}
    If you are logging in for the first time, Seqera sends an authentication link to your email address to complete the login process.
    :::

3. Go to **My Resources** and select the **My access tokens** tab.

4. Under **Create New Access Token**, enter a descriptive name for the **Token name** and select a token duration from the **Expiry** drop down.

5. Select **Generate Token**.

6. Copy the token and store it somewhere safe.
   You cannot view it again after this step.

Once you have your token, see [Publishing modules](#publishing-modules) for instructions on how to use it.

## Viewing modules

View the modules you have published in the registry:

1. Open the [Nextflow module registry](https://registry.nextflow.io/) in a browser.

2. Log in and go to **My Resources**.

3. Select the **My Modules** tab.

The page lists all modules published under your namespaces, along with their versions and status.

(publish-modules)=

## Publishing modules

Publishing requires a {ref}`claimed namespace <module-registry-claim>` and an {ref}`access token <module-registry-access-token>`.

Provide your registry token in one of two ways:

**Environment variable:**

```console
$ export NXF_REGISTRY_TOKEN=<token>
```

**Nextflow configuration:**

```groovy
registry {
    apiKey = '${MYORG_TOKEN}'
}
```

Use the `module publish` command to upload a module to the registry:

```console
$ nextflow module publish myorg/my-module
```

The module directory must include the following files:

- `main.nf`: The module script.
- `meta.yml`: The module spec with name, version, description, and license.
- `README.md`: Module documentation.

Provide a `namespace/name` reference (for an already-installed module) or a local directory path containing the module files.

  :::{tip}
  Before publishing, verify your module is ready:

  - Run {ref}`module validate <cli-module-validate>` to check that all required files and fields are present.
  - Run {ref}`module spec <cli-module-spec>` to generate or update `meta.yml` from your `main.nf`.
  - Add `-dry-run` to `module publish` to validate the structure without uploading.

  See {ref}`dev-modules-page` for a complete guide on creating modules.
  :::

See {ref}`cli-module-publish` for the full command reference.

## Registry configuration

By default, Nextflow uses the public registry at `https://registry.nextflow.io`.
Configure alternative or additional registries using the `registry` scope in your Nextflow configuration.

### Use a private registry

Replace the default registry with your organization's private registry:

```groovy
registry {
    url = 'https://registry.myorg.com'
    apiKey = '${MYORG_TOKEN}'
}
```

Nextflow uses this registry for all module operations (search, install, and publish).

### Use multiple registries

Specify multiple registries as a list.
Nextflow queries them in order when searching for or installing a module:

```groovy
registry {
    url = [
        'https://registry.myorg.com',
        'https://private.registry.nextflow.io'
    ]
    apiKey = '${MYORG_TOKEN}'
}
```

In this example, Nextflow searches the private registry first and falls back to the public registry.

:::{note}
The `apiKey` authenticates with the primary (first) registry.
:::

:::{tip}
Override the target registry with the `-registry` flag:

```console
$ export NXF_REGISTRY_TOKEN=<token>
$ nextflow module publish myorg/my-module -registry 'https://private.registry.myorg.com'
```
:::
