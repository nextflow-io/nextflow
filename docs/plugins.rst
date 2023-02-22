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
When specifying plugins via the command line, any plugin declarations in the Nextflow config file are ignored.


Index
=====

Nextflow resolves plugins download location through the `Plugins index <https://github.com/nextflow-io/plugins/>`_.
The index stores for each plugin the available version, the creation date, checksum and the link from where the plugin
file is downloaded.

To add a new plugin to the Index, create a pull request including the request plugin metadata.
The `nf-hello <https://github.com/nextflow-io/nf-hello>`_ repository provides a minimal code example for
the implementation of a Nextflow plugin.

Import operators from plugin
============================

As of version `22.04.x`, Nextflow allows the inclusion of extension operators from Nextflow plugins.

For example::

    include { sqlInsert; fromQuery as selectFromTable } from 'plugin/nf-sqldb'

    def sql = "select * from FOO"
    channel
        .selectFromTable(sql, db: "test", emitColumns:true)
        .sqlInsert(into:"BAR", columns:'id', db:"test")

The above snippet includes the operators ``sqlInsert`` and ``fromQuery`` from the
`nf-sqldb <https://github.com/nextflow-io/nf-sqldb>`_ plugin. The latter will be accessible using
the ``selectFromTable`` alias in the script.

.. note::
    The prefix ``plugin/`` must precede the plugin name in the include ``from`` statement.


Import custom functions from plugin
===================================

In the same way, as of version `22.09.x`, a plugin can export custom functions.

For example, a plugin can export a util function to reverse a String::

     @nextflow.plugin.extension.Function
     String reverseString( String origin ){
          origin.reverse()
     }

And this function can be used by the pipeline::

    include { reverseString } from 'plugin/my-plugin'

    channel.of( reverseString('hi') )

The above snippet includes a function from the plugin and allows the channel to call it directly.

In the same way as operators, functions can be aliased::

    include { reverseString as anotherReverseMethod } from 'plugin/my-plugin'

