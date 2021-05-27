.. _k8s-page:

**********
Kubernetes
**********

`Kubernetes <https://kubernetes.io/>`_ is a cloud-native open-source system for deployment, scaling, and management of
containerized applications.

It provides clustering and file system abstractions that allows the execution of containerised workloads across
different cloud platforms and on-premises installations.

The built-in support for Kubernetes provided by Nextflow streamlines the execution of containerised workflows in
Kubernetes clusters.


Concepts
========

Kubernetes main abstraction is the `pod`. A `pod` defines the (desired) state of one or more containers i.e. required
computing resources, storage, network configuration.

Kubernetes abstracts also the storage provisioning through the definition of one more more persistent volumes that
allow containers to access to the underlying storage systems in a transparent and portable manner.

When using the ``k8s`` executor Nextflow deploys the workflow execution as a Kubernetes pod. This pod orchestrates
the workflow execution and submits a separate pod execution for each job that need to be carried out by the workflow
application.

.. image:: images/nextflow-k8s-min.png


Requirements
============

At least a `Persistent Volume <https://kubernetes.io/docs/concepts/storage/persistent-volumes/#persistent-volumes>`_ with
``ReadWriteMany`` access mode has to be defined in the Kubernetes cluster (check the supported storage systems
at `this link <https://kubernetes.io/docs/concepts/storage/persistent-volumes/#access-modes>`_).

Such volume needs to be accessible through a
`Persistent Volume Claim <https://kubernetes.io/docs/concepts/storage/persistent-volumes/#persistentvolumeclaims>`_, which
will be used by Nextflow to run the application and store the scratch data and the pipeline final result.

The workflow application has to be containerised using the usual Nextflow :ref:`container<process-container>` directive.


Execution
=========

The workflow execution needs to be submitted from a computer able to connect to the Kubernetes cluster.

Nextflow uses the Kubernetes configuration file available at the path ``$HOME/.kube/config`` or the file specified
by the environment variable ``KUBECONFIG``.

You can verify such configuration with the command below::

    $ kubectl cluster-info
    Kubernetes master is running at https://your-host:6443
    KubeDNS is running at https://your-host:6443/api/v1/namespaces/kube-system/services/kube-dns:dns/proxy


To deploy and launch the workflow execution use the Nextflow command ``kuberun`` as shown below::

    nextflow kuberun <pipeline-name> -v vol-claim:/mount/path


This command will create and execute a pod running the nextflow orchestrator for the specified workflow.
In the above example replace ``<pipeline-name>`` with an existing nextflow project or the absolute path
of a workflow already deployed in the Kubernetes cluster.

The ``-v`` command line option is required to specify the volume claim name and mount path to use for the workflow
execution. In the above example replace ``vol-claim`` with the name of an existing persistent volume claim and
``/mount/path`` with the path where the volume is required to be mount in the container. Volume claims can also be
specified in the Nextflow configuration file, see the :ref:`Kubernetes configuration section<config-k8s>` for details.

Once the pod execution starts, the application in the foreground prints the console output produced by the running
workflow pod.

Interactive login
=================

For debugging purpose it's possible to execute a Nextflow pod and launch an interactive shell using the following command::

   nextflow kuberun login -v vol-claim:/mount/path

This command creates a pod, sets up the volume claim(s), configures the Nextflow environment and finally launches a Bash
login session.  

.. warning:: The pod is automatically destroyed once the shell session terminates. Do not use to start long running
  workflow executions in background.


Running in a pod
==================

The main convenience of the ``kuberun`` command is that it spares the user from manually creating a pod from
where the main Nextflow application is launched. In this scenario, the user environment is not containerised.

However there are scenarios in which Nextflow needs to be executed directly from a pod running in a
Kubernetes cluster. In these cases you will need to use the plain Nextflow ``run`` command and specify
the ``k8s`` executor and the required persistent volume claim in the ``nextflow.config`` file as shown below::

    process {
       executor = 'k8s'
    }

    k8s {
       storageClaimName = 'vol-claim'
       storageMountPath = '/mount/path'
       storageSubPath = '/my-data'
    }

In the above snippet replace ``vol-claim`` with the name of an existing persistent volume claim and replace
``/mount/path`` with the actual desired mount path (default: ``/workspace``) and ``storageSubPath``
with the directory in the volume to be mounted (default: ``/``).

.. warning:: The running pod must have been created with the same persistent volume claim name and mount as the
    one specified in your Nextflow configuration file.
    Note also that the ``run`` command does not support the ``-v`` option.
   

Pod settings
============

The process :ref:`process-pod` directive allows the definition of pods specific settings, such as environment variables,
secrets and config maps when using the :ref:`k8s-executor` executor. See the :ref:`process-pod` directive for more details.

Limitation
==========

Currently, the ``kuberun`` command does not allow the execution of local Nextflow scripts.


Advanced configuration
======================

Read :ref:`Kubernetes configuration<config-k8s>` and :ref:`executor <k8s-executor>` sections to learn more
about advanced configuration options.
