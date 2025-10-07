(plugin-registry-page)=

# Nextflow plugin registry

The [Nextflow plugin registry](https://registry.nextflow.io/) is a central repository for Nextflow plugins. It hosts an index of plugin metadata that supports plugin discovery, accessibility, and version tracking.

Nextflow 25.10 and later can use the plugin registry as a drop-in replacement for the [legacy plugin index](https://github.com/nextflow-io/plugins) hosted on GitHub. See {ref}`migrate-plugin-page` for more information about migrating to the Nextflow plugin registry.

:::{note}
The Nextflow plugin registry is currently available as a public preview. Plugin developers can access the registry by contacting [info@nextflow.io](mailto:info@nextflow.io) to obtain access to the registry.
:::

(plugin-registry-claim)=

## Claiming a plugin

Ownership of a plugin is required to publish plugins to the Nextflow plugin registry.

To claim ownership of a plugin:

1. Open the [Nextflow plugin registry](https://registry.nextflow.io/) in a browser.

2. Log in to [Seqera](https://cloud.seqera.io/login) with your GitHub or Google account, or by providing an email address.

    :::{note}
    If you are logging in for the first time, Seqera will send an authentication link to your email address to complete the login process.
    :::

3. Go to the **My plugins** page and select **Claim a plugin**.

4. Enter your unique plugin name or select the plugin you wish to claim in the **Plugin name** field.

5. Enter your organization name in the **Provider** field.

    :::{note}
    This organization must match with the organization specified when publishing your plugin.
    :::

6. Select **Submit Request**.

The plugin will show as **PENDING REVIEW** under **Pending Ownership Requests** until an admin approves the claim. Admin approval is required only once.

(plugin-registry-access-token)=

## Creating an access token

An API access token is required to publish plugins to the Nextflow plugin registry.

To create an API access token:

1. Open the [Nextflow plugin registry](https://registry.nextflow.io/) in a browser.

2. Log in to [Seqera](https://cloud.seqera.io/login) with your GitHub or Google account, or by providing an email address.

    :::{note}
    If you are logging in for the first time, Seqera will send an authentication link to your email address to complete the login process.
    :::

3. Go to the **Access tokens** page.

4. Under **Create New Access Token**, enter a descriptive name for the **Token name** and select the token duration from the **Expiry** drop down.

5. Select **Generate token**.

6. Copy and past token somewhere safe, you won't be able to see it again.

Once you have your token, see {ref}`gradle-plugin-publish` for instructions on how to use it.
