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

Nextflow *defaults* plugins do not require special configuration, they are automatically installed on-demand when
the corresponding feature is requested by a Nextflow pipeline.

Non-defaults plugins needs to be declared in the ``nextflow.config`` file to make them available during the pipeline execution
as shown below:

```
plugins {
  id 'nf-ignite@1.1.0'
}
```



Registry
=========

Plugins need to be listed in the [Plugins registry](https://github.com/nextflow-io/plugins/blob/main/plugins.json) to be
accessible from Nextflow. The registry stores for each plugin the available version, the creation date, checksum and the
link from where the plugin file con be download.

