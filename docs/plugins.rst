.. _plugins-page:

*********
Plugins
*********

Main concepts
=============

Nextflow is based on a plugins system that allows extending core functionalities via pluggable components
that are download and installed at runtime.

Currently the following functionalities are implemented as plugins component and makes part of the
Nextflow default plugins:

* ``nf-amazon``: Support for Amazon cloud.
* ``nf-azure``: Support for Azure cloud.
* ``nf-console``: Implement Nextflow `REPL console <https://www.nextflow.io/blog/2015/introducing-nextflow-console.html>`_.
* ``nf-ga4gh``: Support `GA4GH APIs <https://www.ga4gh.org/>`_.
* ``nf-google``: Support for Google cloud.
* ``nf-ignite``: Implement Nextflow :ref:`Ignite executor<ignite-page>`.
* ``nf-tower``: Support for `Nextflow Tower <https://tower.nf>`_ platform.

Nextflow defaults plugins do not require special configuration, they are automatically installed on-demand when
the corresponding feature is use by a Nextflow pipeline.

Registry
=========
