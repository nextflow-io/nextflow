.. _fusion-page:

******************
Fusion file system
******************

Introduction
=============

Fusion is a distributed virtual file system for cloud-native data pipeline and optimised for Nextflow workloads.

It bridges the gap between cloud-native storage and data analysis workflow by implementing a thin client
that allows any existing application to access object storage using the standard POSIX interface, thus simplifying
and speeding up most operations. Currently it supports AWS S3.

.. warning:: This is an incubating feature. It may change in future Nextflow releases.

Getting started
===============

Requirements
-------------

Fusion file system is designed to work with containerised workloads, therefore it requires the use of a container
engine such as Docker or a container native platform for the execution of your pipeline e.g. AWS Batch or Kubernetes.

It also requires the use of :ref:`Wave containers<wave-page>` and Nextflow version ``22.10.0`` or later.

AWS S3 configuration
--------------------

The AWS S3 bucket should be configured with the following IAM permissions::

    {
        "Version": "2012-10-17",
        "Statement": [
            {
                "Effect": "Allow",
                "Action": [
                    "s3:ListBucket"
                ],
                "Resource": [
                    "arn:aws:s3:::YOUR-BUCKET-NAME"
                ]
            },
            {
                "Action": [
                    "s3:GetObject",
                    "s3:PutObject",
                    "s3:PutObjectTagging",
                    "s3:DeleteObject"
                ],
                "Resource": [
                    "arn:aws:s3:::YOUR-BUCKET-NAME/*"
                ],
                "Effect": "Allow"
            }
        ]
    }


Use cases
=========

Local execution with S3 bucket as work directory
------------------------------------------------

Fusion file system allows the use of an S3 bucket as a pipeline work directory with the Nextflow local executor. This
configuration requires the use of Docker (or similar container engine) for the execution of your pipeline tasks.

The AWS S3 bucket credentials should be made accessible via standard ``AWS_ACCESS_KEY_ID`` and ``AWS_SECRET_ACCESS_KEY``
environment variables.

The following configuration should be added in your Nextflow configuration file::

    docker {
      enabled = true
    }

    fusion {
      enabled = true
      exportAwsAccessKeys = true
    }

    wave {
        enabled = true
    }


Then you can run your pipeline using the following command::

    nextflow run <YOUR PIPELINE> -work-dir s3://<YOUR BUCKET>/scratch

Replace ``<YOUR PIPELINE>`` and ``<YOUR BUCKET>`` with a pipeline script and bucket or your choice, for example::

    nextflow run https://github.com/nextflow-io/rnaseq-nf -work-dir s3://nextflow-ci/scratch


AWS Batch execution with S3 bucket as work directory
----------------------------------------------------

Fusion file system allows the use of an S3 bucket as a pipeline work directory with the AWS Batch executor. The use
of Fusion makes obsolete the need to create and configure a custom AMI that includes the `aws` command line tool, when
setting up the AWS Batch compute environment.

The configuration for this deployment scenario looks like the following::

    fusion {
      enabled = true
    }

    wave {
      enabled = true
    }

    process {
      executor = 'awsbatch'
      queue = '<YOUR BATCH QUEUE>'
    }

    aws {
      region = '<YOUR AWS REGION>'
    }

Then you can run your pipeline using the following command::

    nextflow run <YOUR PIPELINE> -work-dir s3://<YOUR BUCKET>/scratch



Kubernetes execution with S3 bucket as work directory
-----------------------------------------------------

Fusion file system allows the use of an S3 bucket as a pipeline work directory with the Kubernetes executor.

The use of Fusion makes obsolete the need to create and manage and separate persistent volume and shared file system
in the Kubernetes cluster.

The configuration for this deployment scenario looks like the following::

    wave {
      enabled = true
    }

    fusion {
      enabled = true
    }

    process {
      executor = 'k8s'
    }

    k8s {
      context = '<YOUR K8S CONFIGURATION CONTEXT>'
      namespace = '<YOUR K8S NAMESPACE>'
      serviceAccount = '<YOUR K8S SERVICE ACCOUNT>'
    }


The ``k8s.context`` represents the Kubernetes configuration context to be used for the pipeline execution. This
setting can be omitted if Nextflow itself is run as a pod in the Kubernetes clusters.

The ``k8s.namespace`` represents the Kubernetes namespace where the jobs submitted by the pipeline execution should
be executed.

The ``k8s.serviceAccount`` represents the Kubernetes service account that should be used to grant the execution
permission to jobs launched by Nextflow. You can find more details how to configure it as the `following link <https://github.com/seqeralabs/wave-showcase/tree/master/example8>`_.


Having the above configuration in place, you can run your pipeline using the following command::

    nextflow run <YOUR PIPELINE> -work-dir s3://<YOUR BUCKET>/scratch

More examples
=============

Check out the `Wave showcase repository <https://github.com/seqeralabs/wave-showcase>`_ for more examples on how to use
Fusion file system.
