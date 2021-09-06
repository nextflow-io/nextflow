.. _amazons3-page:

*******************
Amazon S3 storage
*******************

Nextflow includes the support for Amazon S3 storage. Files stored in a S3 bucket can be accessed
transparently in your pipeline script like any other file in the local file system.

S3 path
---------
In order to access a S3 file you only need to prefix the file path with the ``s3`` schema and the `bucket` name
where it is stored.

For example if you need to access the file ``/data/sequences.fa`` stored in a bucket with name ``my-bucket``,
that file can be accessed using the following fully qualified path::

   s3://my-bucket/data/sequences.fa


The usual file operations can be applied on a path handle created using the above notation. For example the content
of a S3 file can be printed as shown below::

    println file('s3://my-bucket/data/sequences.fa').text


See section :ref:`script-file-io` to learn more about available file operations.




Security credentials
---------------------

Amazon access credentials can be provided in two ways:

#. Using AWS access and secret keys in your pipeline configuration.
#. Using IAM roles to grant access to S3 storage on Amazon EC2 instances.

AWS access and secret keys
===========================

The AWS access and secret keys can be specified by using the ``aws`` section in the ``nextflow.config`` configuration
file as shown below::

  aws {
    accessKey = '<Your AWS access key>'
    secretKey = '<Your AWS secret key>'
    region = '<AWS region identifier>'
  }


If the access credentials are not found in the above file, Nextflow looks for AWS credentials in a number of different
places, including environment variables and local AWS configuration files.


Nextflow looks for AWS credentials in the following order:

    #. the ``nextflow.config`` file in the pipeline execution directory
    #. the environment variables ``AWS_ACCESS_KEY_ID`` and ``AWS_SECRET_ACCESS_KEY``
    #. the environment variables ``AWS_ACCESS_KEY`` and ``AWS_SECRET_KEY``
    #. the `default` profile in the AWS credentials file located at ``~/.aws/credentials``
    #. the `default` profile in the AWS client configuration file located at ``~/.aws/config``
    #. the temporary AWS credentials provided by an IAM instance role. See `IAM Roles <http://docs.aws.amazon.com/AWSEC2/latest/UserGuide/iam-roles-for-amazon-ec2.html>`_ documentation for details.


More information regarding `AWS Security Credentials <http://docs.aws.amazon.com/general/latest/gr/aws-security-credentials.html>`_
are available in Amazon documentation.

IAM roles Amazon EC2 instances
================================

When running your pipeline into a Ec2 instance, IAM roles can be used to grant access to AWS resources.

In this scenario, you only need to launch the Ec2 instance specifying a IAM role which includes a
`S3 full access` policy. Nextflow will detected and acquire automatically the access grant to the S3 storage,
without any further configuration.

Learn more about `Using IAM Roles to Delegate Permissions to Applications that Run on Amazon EC2 <http://docs.aws.amazon.com/IAM/latest/UserGuide/roles-usingrole-ec2instance.html>`_ on Amazon
documentation.

Advanced configuration
-----------------------

Read :ref:`AWS configuration<config-aws>` section to learn more about advanced S3 client configuration options.







