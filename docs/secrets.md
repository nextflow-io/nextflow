(secrets-page)=

# Secrets

:::{versionadded} 22.10.0
:::

Nextflow has built-in support for pipeline secrets, allowing users to safely provide sensitive information to a pipeline run.

## How it works

This feature allows you to decouple the use of secrets in your pipelines from the pipeline code and configuration files. Secrets are managed by Nextflow and stored separately into a local store only accessible to the secrets owner.

When a pipeline is launched, Nextflow injects the secrets into the run without leaking them into temporary execution files. Secrets are provided to tasks as environment variables.

Secrets can be used with the local executor and grid executors (e.g., Slurm or Grid Engine). Secrets can be used with the AWS Batch executor when launched from [Seqera Platform](https://seqera.io/blog/pipeline-secrets-secure-handling-of-sensitive-information-in-tower/).

## Command line

The Nextflow {ref}`cli-secrets` sub-command can be used to manage secrets:

```bash
# create a new secret
nextflow secrets set FOO "Hello world"

# list all secrets
nextflow secrets list

# get the value of secret FOO
nextflow secrets get FOO

# delete the secret FOO
nextflow secrets delete FOO
```

## Configuration file

Secrets can be used in configuration files using the built-in `secrets` variable. For example:

```groovy
aws {
  accessKey = secrets.MY_ACCESS_KEY
  secretKey = secrets.MY_SECRET_KEY
}
```

The above snippet accesses the secrets `MY_ACCESS_KEY` and `MY_SECRET_KEY` and assigns them to the corresponding AWS config settings.

:::{warning}
Secrets should not be assigned to pipeline parameters, as they can be leaked by the pipeline.
:::

:::{versionadded} 25.10.0
:::

Nextflow supports the use of secrets provided by plugins (e.g., AWS secrets) in configuration. However, due to the way that plugins are loaded, there are specific considerations when using config secrets:

- **Initial config load**: Nextflow first loads the configuration _without_ secrets enabled. Any reference to a secret will return the empty string `''`.

- **Plugin resolution**: Plugins are resolved after the initial configuration load. This is because the configuration can specify additional plugins.

- **Config reloading**: If secrets are accessed during configuration and the initial load succeeds, Nextflow will reload the configuration with secrets enabled.

As a result, config secrets must be used in a way that does not cause the config resolution to fail when secrets are not present.

For example:

```groovy
includeConfig secrets.MY_SECRET
    ? "https://example.com/extra.config?secret=${secrets.MY_SECRET}"
    : '/dev/null'
```

The above snippet includes a secured config only if the secret is present. Otherwise, it includes `/dev/null`, which is equivalent to including an empty file. The reference to `secrets.MY_SECRET` in the condition causes the config to be reloaded with secrets enabled, including secrets from plugins such as AWS secrets.

(secrets-pipeline-script)=

## Pipeline script

:::{versionadded} 24.04.0
:::

Secrets can be accessed in the pipeline script using the built-in `secrets` variable. For example:

```nextflow
workflow {
    println "The secret is: ${secrets.MY_SECRET}"
}
```

:::{warning}
The above example is only meant to demonstrate how to access a secret, not how to use it. In practice, sensitive information should not be printed to the console or output files.
:::

## Process directive

Secrets can be accesses by processes using the {ref}`process-secret` directive. For example:

```nextflow
process my_task {
    secret 'MY_ACCESS_KEY'
    secret 'MY_SECRET_KEY'

    script:
    """
    your_command --access \$MY_ACCESS_KEY --secret \$MY_SECRET_KEY
    """
}
```

In the above example, the secrets `MY_ACCESS_KEY` and `MY_SECRET_KEY` are injected into the process script as environment variables.

:::{warning}
Secrets are made available as environment variables in the process script. To prevent evaluation in the Nextflow script context, escape variable names with a backslash (e.g., `\$MY_ACCESS_KEY`) as shown above.
:::
