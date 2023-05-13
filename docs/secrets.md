(secrets-page)=

# Secrets

:::{versionadded} 22.10.0
Previewed in `21.09.0-edge`.
:::

Nextflow has built-in support for pipeline secrets to allow users to safely provide sensitive information to a pipeline execution.

## How it works

This feature allows decoupling the use secrets in your pipelines from the pipeline code and configuration files. Secrets are instead managed by Nextflow and store separately into a local store only accessible to the secrets owner.

When the pipeline execution is launched Nextflow inject the secrets in pipeline jobs without leaking them into temporary execution files. The secrets are accessible into the job command via environment variables.

## Command line

Nextflow provides a command named `secrets`. This command allows four simple operations:

`list`
: List secrets available in the current store e.g. `nextflow secrets list`.

`get`
: Retrieve a secret value e.g. `nextflow secrets get FOO`.

`set`
: Create or update a secret e.g. `nextflow secrets set FOO "Hello world"`

`delete`
: Delete a secret e.g. `nextflow secrets delete FOO`.

## Configuration file

Once create the secrets can be used in the pipeline configuration file as implicit variables using the `secrets` scope:

```groovy
aws {
  accessKey = secrets.MY_ACCESS_KEY
  secretKey = secrets.MY_SECRET_KEY
}
```

The above snippet access the secrets `MY_ACCESS_KEY` and `MY_SECRET_KEY` previously and assign them to the corresponding AWS credentials settings.

:::{warning}
Secrets **cannot** be assigned to pipeline parameters.
:::

## Process secrets

Secrets can be access by pipeline processes by using the `secret` directive. For example:

```groovy
process someJob {
    secret 'MY_ACCESS_KEY'
    secret 'MY_SECRET_KEY'

    """
    your_command --access \$MY_ACCESS_KEY --secret \$MY_SECRET_KEY
    """
}
```

The above snippet runs a command in with the variables `MY_ACCESS_KEY` and `MY_SECRET_KEY` are injected in the process execution environment holding the values defines in the secret store.

:::{warning}
The secrets are made available in the process context running the command script as environment variables. Therefore make sure to escape the variable name identifier with a backslash as shown in the example above, otherwise a variable with the same will be evaluated in the Nextflow script context instead of the command script.
:::

:::{note}
This feature is only available when using the local or grid executors (Slurm, Grid Engine, etc). The AWS Batch executor allows the use of secrets when deploying the pipeline execution via [Nextflow Tower](https://seqera.io/blog/pipeline-secrets-secure-handling-of-sensitive-information-in-tower/).
:::
