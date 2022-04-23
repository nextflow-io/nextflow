.. _plugins-page:

*******
Plugins
*******

Main concepts
=============

Nextflow is based on a plugins system that allows extending core functionalities via pluggable components
that are download and installed at runtime.

Currently the following functionalities are implemented as plugin components and they make part of the
Nextflow *default* plugins:

* ``nf-amazon``: Support for Amazon cloud.
* ``nf-azure``: Support for Azure cloud.
* ``nf-console``: Implement Nextflow `REPL console <https://www.nextflow.io/blog/2015/introducing-nextflow-console.html>`_.
* ``nf-ga4gh``: Support `GA4GH APIs <https://www.ga4gh.org/>`_.
* ``nf-google``: Support for Google cloud.
* ``nf-tower``: Support for `Nextflow Tower <https://tower.nf>`_ platform.


Configuration
==============

Nextflow *default* plugins do not require any configuration. They are automatically installed when
the corresponding feature is requested by a Nextflow pipeline.

To use **non-default** plugins in your pipeline execution, you must declare them in the Nextflow configuration file,
listing each plugin as shown below::

    plugins {
      id 'nf-hello@0.1.0'
    }

The plugin identifier consists of the plugin name and plugin version separated by a ``@``.

Alternatively, plugins can be required using the ``-plugins`` command line option::

    nextflow run <PIPELINE NAME> -plugins nf-hello@0.1.0

Multiple plugins can be specified by separating them with a comma.
When specifiying plugins via the command line, any plugin declarations in the Nextflow config file are ignored.


Index
=====

Nextflow resolves plugins download location through the `Plugins index <https://github.com/nextflow-io/plugins/>`_.
The index stores for each plugin the available version, the creation date, checksum and the link from where the plugin
file is downloaded.

To add a new plugin to the Index, create a pull request including the request plugin metadata.
The `nf-hello <https://github.com/nextflow-io/nf-hello>`_ repository provides a minimal code example for
the implementation of a Nextflow plugin.
