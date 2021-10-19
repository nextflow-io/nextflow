.. _secrets-page:

*******
Secrets
*******


As of version ``21.09.0-edge``, Nextflow adds the built-in support for pipeline secrets to allow users to handle
and manage sensitive information for pipeline execution in a safe manner.

.. warning::
    This is a preview feature, therefore options and syntax may change in future release.

How it works
============

This feature allows decoupling the use secrets in your pipelines from the pipeline code and configuration files.
Secrets are instead managed by Nextflow and store separately into a local store only accessible to the secrets
owner.

When the pipeline execution is launched Nextflow inject the secrets in pipeline jobs without leaking them
into temporary execution files. The secrets are accessible into the job command via environment variables.

.. note::
  This feature needs to be enabled by settings the following environment variable in the launching environment::

        export NXF_ENABLE_SECRETS=true


Command line
============

When enabling this feature Nextflow provides a new command named ``secrets``. This command allow four simple
operations:

===================== =====================
Operation               Description
===================== =====================
``list``                List secrets available in the current store e.g. ``nextflow secrets list``.
``get``                 Allows retrieving a secret value e.g. ``nextflow secrets get -n FOO``.
``put``                 Allows creating creating a new secret or overriding an existing one e.g. ``nextflow secrets put -n FOO -v "Hello world"``
``delete``              Allows deleting an existing secret e.g. ``nextflow secrets delete -n FOO``.
===================== =====================

Configuration file
==================

Once create the secrets can be used in the pipeline configuration file as implicit variables using the ``secrets`` scope::

    aws {
      accessKey = secrets.MY_ACCESS_KEY
      secretKey = secrets.MY_SECRET_KEY
    }

The above above snippet access the secrets ``MY_ACCESS_KEY`` and ``MY_SECRET_KEY`` previously and assign them to
the corresponding AWS credentials settings.

.. note::
    Secrets **cannot** be assigned to pipeline parameters. 


Process secrets
===============

Secrets can be access by pipeline processes by using the `secret` directive. For example::

    process someJob {
        secret 'MY_ACCESS_KEY'
        secret 'MY_SECRET_KEY'

        """
        your_command --access $MY_ACCESS_KEY --secret $MY_SECRET_KEY
        """
    }

The above snippet runs a command in with the variables ``MY_ACCESS_KEY`` and ``MY_SECRET_KEY`` are injected in the
process execution environment holding the values defines in the secret store.

.. warning::
    This feature is only available when using the local or batch scheduler executions e.g. Slurm, Grid Engine, etc.
