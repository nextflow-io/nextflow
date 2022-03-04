.. _ga4ghtes-executor:

GA4GH TES
=========

.. warning:: This is an experimental feature and it may change in a future release. It requires Nextflow
  version 0.31.0 or later.

The `Task Execution Schema <https://github.com/ga4gh/task-execution-schemas>`_ (TES) project
by the `GA4GH <https://www.ga4gh.org>`_ standardisation initiative is an effort to define a
standardized schema and API for describing batch execution tasks in portable manner.

Nextflow includes an experimental support for the TES API providing a ``tes`` executor which allows
the submission of workflow tasks to a remote execution back-end exposing a TES API endpoint.

To use this feature define the following variables in the workflow launching environment::

    export NXF_MODE=ga4gh
    export NXF_EXECUTOR=tes
    export NXF_EXECUTOR_TES_ENDPOINT='http://back.end.com'

It is important that the endpoint is specified without the trailing slash; otherwise, the resulting URLs will be not
normalized and the requests to TES will fail.

Then you will be able to run your workflow over TES using the usual Nextflow command line. Be sure to specify the Docker
image to use, i.e.::

    nextflow run rnaseq-nf -with-docker alpine

.. note:: If the variable ``NXF_EXECUTOR_TES_ENDPOINT`` is omitted the default endpoint is ``http://localhost:8000``.

.. tip:: You can use a local `Funnel <https://ohsu-comp-bio.github.io/funnel/>`_ server using the following launch
  command line::

  ./funnel server --Server.HTTPPort 8000 --LocalStorage.AllowedDirs $HOME run

  (tested with version 0.8.0 on macOS)

.. warning:: Make sure the TES back-end can access the workflow work directory when
  data is exchanged using a local or shared file system.

**Known Limitations**

* Automatic deployment of workflow scripts in the `bin` folder is not supported.
* Process output directories are not supported. For details see `#76 <https://github.com/ga4gh/task-execution-schemas/issues/76>`_.
* Glob patterns in process output declarations are not supported. For details see `#77 <https://github.com/ga4gh/task-execution-schemas/issues/77>`_.
