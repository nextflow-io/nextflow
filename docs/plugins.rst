.. _plugins-page:

*********
Plugins
*********

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

Nextflow *defaults* plugins do not require any configuration, they are automatically installed on-demand when
the corresponding feature is requested by a Nextflow pipeline.

To use **non-default** plugins in your pipeline execution it's required to declared them in the Nextflow configuration file
listing each plugin identifiers as shown below::

    plugins {
      id 'nf-hello@0.1.0'
    }


.. note::
  The plugin identifier is composed by the plugin name, followed by the ``@`` separator and finally the plugin version.

Alternatively, plugins can be required using the command line option ``-plugins`` e.g.::

    nextflow run <PIPELINE NAME> -plugins nf-hello@0.1.0


.. note::
  When using the ``-plugins`` CLI option any plugin declaration in the Nextflow config file is ignored.
  Multiple plugin Ids can be specified separating them with a comma character.


Index
======

Nextflow resolves plugins download location through the `Plugins index <https://github.com/nextflow-io/plugins/>`_.
The index stores for each plugin the available version, the creation date, checksum and the link from where the plugin
file is downloaded.

To add a new plugin to the Index, create a pull request including the request plugin metadata.

.. tip::
  The `nf-hello plugin <https://github.com/nextflow-io/nf-hello>`_ repository provides an bare minimal code example for
  the implementation of a Nextflow plugin.
