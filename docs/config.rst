.. _config-page:

*************
Configuration
*************

Configuration file
==================

When a pipeline script is launched, Nextflow looks for configuration files in multiple locations.
Since each configuration file can contain conflicting settings, the sources are ranked to decide which
settings to are applied. All possible configuration sources are reported below, listed in order
of priority:

1. Parameters specified on the command line (``--something value``)
2. Parameters provided using the ``-params-file`` option
3. Config file specified using the ``-c my_config`` option
4. The config file named ``nextflow.config`` in the current directory
5. The config file named ``nextflow.config`` in the workflow project directory
6. The config file ``$HOME/.nextflow/config``
7. Values defined within the pipeline script itself (e.g. ``main.nf``)

When more than one of these ways of specifying configurations are used, they are merged, so that the settings in the
first override the same ones that may appear in the second one, and so on.

.. tip::
  If you want to ignore any default configuration files and use only the custom one, use ``-C <config file>``.


Config syntax
-------------

A Nextflow configuration file is a simple text file containing a set of properties defined using the syntax::

  name = value

Please note, string values need to be wrapped in quotation characters while numbers and boolean values (``true``, ``false``) do not.
Also note that values are typed, meaning for example that, ``1`` is different from ``'1'``, since the first is interpreted
as the number one, while the latter is interpreted as a string value.


Config variables
----------------

Configuration properties can be used as variables in the configuration file itself, by using the usual
``$propertyName`` or ``${expression}`` syntax.

For example::

     propertyOne = 'world'
     anotherProp = "Hello $propertyOne"
     customPath = "$PATH:/my/app/folder"

Please note, the usual rules for :ref:`string-interpolation` are applied, thus a string containing a variable
reference must be wrapped in double-quote chars instead of single-quote chars.

The same mechanism allows you to access environment variables defined in the hosting system. Any variable whose name is
not defined in the Nextflow configuration file(s) is supposed to be a reference to an environment variable with that name.
So, in the above example the property ``customPath`` is defined as the current system ``PATH`` to which
the string ``/my/app/folder`` is appended.


Config comments
---------------

Configuration files use the same conventions for comments used by the Groovy or Java programming languages. Thus, use ``//`` to comment
a single line or ``/*`` .. ``*/`` to comment a block on multiple lines.


Config include
--------------

A configuration file can include one or more configuration files using the keyword ``includeConfig``. For example::

    process.executor = 'sge'
    process.queue = 'long'
    process.memory = '10G'

    includeConfig 'path/foo.config'

When a relative path is used, it is resolved against the actual location of the including file.


Config scopes
=============

Configuration settings can be organized in different scopes by dot prefixing the property names with a scope
identifier or grouping the properties in the same scope using the curly brackets notation. This is shown in the
following example::

   alpha.x  = 1
   alpha.y  = 'string value..'

   beta {
        p = 2
        q = 'another string ..'
   }


.. _config-aws:

Scope `aws`
-----------

The ``aws`` scope allows you to configure the access to Amazon S3 storage. Use the attributes ``accessKey`` and ``secretKey``
to specify your bucket credentials. For example::


    aws {
        accessKey = '<YOUR S3 ACCESS KEY>'
        secretKey = '<YOUR S3 SECRET KEY>'
        region = '<REGION IDENTIFIER>'
    }

Click the following link to learn more about `AWS Security Credentials <http://docs.aws.amazon.com/general/latest/gr/aws-security-credentials.html>`_.

Advanced client configuration options can be set by using the ``client`` attribute. The following properties can be used:

=========================== ================
Name                        Description
=========================== ================
anonymous                   Allow the access of public S3 buckets without the need to provide AWS credentials. Any service that does not accept unsigned requests will return a service access error.
s3Acl                       Allow the setting of a predefined bucket permissions also known as *canned ACL*. Permitted values are ``Private``, ``PublicRead``, ``PublicReadWrite``, ``AuthenticatedRead``, ``LogDeliveryWrite``, ``BucketOwnerRead``, ``BucketOwnerFullControl`` and ``AwsExecRead``. See `Amazon docs <https://docs.aws.amazon.com/AmazonS3/latest/userguide/acl-overview.html#canned-acl>`_ for details.
connectionTimeout           The amount of time to wait (in milliseconds) when initially establishing a connection before giving up and timing out.
endpoint                    The AWS S3 API entry point e.g. `s3-us-west-1.amazonaws.com`.
maxConnections              The maximum number of allowed open HTTP connections.
maxErrorRetry               The maximum number of retry attempts for failed retryable requests.
protocol                    The protocol (i.e. HTTP or HTTPS) to use when connecting to AWS.
proxyHost                   The proxy host to connect through.
proxyPort                   The port on the proxy host to connect through.
proxyUsername               The user name to use when connecting through a proxy.
proxyPassword               The password to use when connecting through a proxy.
s3PathStyleAccess           Enable the use of path-based access model that is used to specify the address of an object in S3-compatible storage systems.
signerOverride              The name of the signature algorithm to use for signing requests made by the client.
socketSendBufferSizeHint    The Size hint (in bytes) for the low level TCP send buffer.
socketRecvBufferSizeHint    The Size hint (in bytes) for the low level TCP receive buffer.
socketTimeout               The amount of time to wait (in milliseconds) for data to be transferred over an established, open connection before the connection is timed out.
storageEncryption           The S3 server side encryption to be used when saving objects on S3, either ``AES256`` or ``aws:kms`` values are allowed.
storageKmsKeyId             The AWS KMS key Id to be used to encrypt files stored in the target S3 bucket (requires version ``22.05.0-edge`` or later).
userAgent                   The HTTP user agent header passed with all HTTP requests.
uploadMaxThreads            The maximum number of threads used for multipart upload.
uploadChunkSize             The size of a single part in a multipart upload (default: `100 MB`).
uploadStorageClass          The S3 storage class applied to stored objects, one of [`STANDARD`, `STANDARD_IA`, `ONEZONE_IA`, `INTELLIGENT_TIERING`] (default: `STANDARD`).
uploadMaxAttempts           The maximum number of upload attempts after which a multipart upload returns an error (default: `5`).
uploadRetrySleep            The time to wait after a failed upload attempt to retry the part upload (default: `500ms`).
=========================== ================

For example::

    aws {
        client {
            maxConnections = 20
            connectionTimeout = 10000
            uploadStorageClass = 'INTELLIGENT_TIERING'
            storageEncryption = 'AES256'
        }
    }


.. _config-aws-batch:

Advanced Batch configuration options can be set by using the ``batch`` attribute. The following properties can be used (required version `19.07.0` or later):

=========================== ================
Name                        Description
=========================== ================
cliPath                     The path where the AWS command line tool is installed in the host AMI.
jobRole                     The AWS Job Role ARN that needs to be used to execute the Batch Job.
logsGroup                   The name of the logs group used by Batch Jobs (default: ``/aws/batch``, requires ``22.09.0-edge`` or later).
volumes                     One or more container mounts. Mounts can be specified as simple e.g. `/some/path` or canonical format e.g. ``/host/path:/mount/path[:ro|rw]``. Multiple mounts can be specifid separating them with a comma or using a list object.
delayBetweenAttempts        Delay between download attempts from S3 (default `10 sec`).
maxParallelTransfers        Max parallel upload/download transfer operations *per job* (default: ``4``).
maxTransferAttempts         Max number of downloads attempts from S3 (default: `1`).
maxSpotAttempts             Max number of execution attempts of a job interrupted by a EC2 spot reclaim event (default: ``5``, requires ``22.04.0`` or later)
shareIdentifier             The share identifier for all tasks when using `fair-share scheduling for AWS Batch <https://aws.amazon.com/blogs/hpc/introducing-fair-share-scheduling-for-aws-batch/>`_ (requires ``22.09.0-edge`` or later)
retryMode                   The retry mode configuration setting, to accommodate rate-limiting on `AWS services <https://docs.aws.amazon.com/cli/latest/userguide/cli-configure-retries.html>`_ (default: ``standard``)
=========================== ================


.. _config-charliecloud:

Scope `charliecloud`
--------------------

The ``charliecloud`` configuration scope controls how `Charliecloud <https://hpc.github.io/charliecloud/>`_ containers are executed by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Charliecloud execution (default: ``false``).
envWhitelist        Comma separated list of environment variable names to be included in the container environment.
temp                Mounts a path of your choice as the ``/tmp`` directory in the container. Use the special value ``auto`` to create a temporary directory each time a container is created.
runOptions          This attribute can be used to provide any extra command line options supported by the ``ch-run`` command.
cacheDir            The directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.
pullTimeout         The amount of time the Charliecloud pull can last, exceeding which the process is terminated (default: ``20 min``).
================== ================

The above options can be used by prefixing them with the ``charliecloud`` scope or surrounding them by curly
brackets, as shown below::

    process.container = 'nextflow/examples'

    charliecloud {
        enabled = true
    }

Read :ref:`container-charliecloud` page to learn more about how to use Charliecloud containers with Nextflow.


.. _config-cloud:

Scope `cloud`
-------------

.. note::
    The ``cloud`` configuration scope is no longer used. See the platform-specific cloud executors instead.


.. _config-conda:

Scope `conda`
-------------

The ``conda`` scope allows for the definition of the configuration settings that control the creation of a Conda environment
by the Conda package manager.

The following settings are available:

================== ================
Name                Description
================== ================
cacheDir            Defines the path where Conda environments are stored. When using a compute cluster make sure to provide a shared file system path accessible from all compute nodes.
createOptions       Defines any extra command line options supported by the ``conda create`` command. For details `Conda documentation <https://docs.conda.io/projects/conda/en/latest/commands/create.html>`_.
createTimeout       Defines the amount of time the Conda environment creation can last. The creation process is terminated when the timeout is exceeded (default: ``20 min``).
useMamba            Uses the ``mamba`` binary instead of ``conda`` to create the Conda environments. For details `Mamba documentation <https://github.com/mamba-org/mamba>`_.
useMicromamba       uses the ``micromamba`` binary instead of ``conda`` to create the Conda environments (requires version ``22.05.0-edge`` or later). For details see `Micromamba documentation <https://mamba.readthedocs.io/en/latest/user_guide/micromamba.html>`_.
================== ================


.. _config-dag:

Scope `dag`
-------------

The ``dag`` scope allows you to control the layout of the execution graph file generated by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             When ``true`` turns on the generation of the execution graph report file (default: ``false``).
file                Graph file name (default: ``dag.dot``).
================== ================

The above options can be used by prefixing them with the ``dag`` scope or surrounding them by curly
brackets. For example::

    dag {
        enabled = true
        file = 'pipeline_dag.html'
    }

To learn more about the execution graph that can be generated by Nextflow read :ref:`dag-visualisation` documentation page.


.. _config-docker:

Scope `docker`
--------------

The ``docker`` configuration scope controls how `Docker <https://www.docker.com>`_ containers are executed by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Docker execution (default: ``false``).
envWhitelist        Comma separated list of environment variable names to be included in the container environment.
legacy              Uses command line options removed since version 1.10.x (default: ``false``).
sudo                Executes Docker run command as ``sudo`` (default: ``false``).
tty                 Allocates a pseudo-tty (default: ``false``).
temp                Mounts a path of your choice as the ``/tmp`` directory in the container. Use the special value ``auto`` to create a temporary directory each time a container is created.
remove              Clean-up the container after the execution (default: ``true``). For details see: https://docs.docker.com/engine/reference/run/#clean-up---rm .
runOptions          This attribute can be used to provide any extra command line options supported by the ``docker run`` command. For details see: https://docs.docker.com/engine/reference/run/ .
registry            The registry from where Docker images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. ``http://``.
fixOwnership        Fixes ownership of files created by the docker container.
engineOptions       This attribute can be used to provide any option supported by the Docker engine i.e. ``docker [OPTIONS]``.
mountFlags          Add the specified flags to the volume mounts e.g. `mountFlags = 'ro,Z'`
================== ================

The above options can be used by prefixing them with the ``docker`` scope or surrounding them by curly
brackets, as shown below::

    process.container = 'nextflow/examples'

    docker {
        enabled = true
        temp = 'auto'
    }

Read :ref:`container-docker` page to learn more about how to use Docker containers with Nextflow.


.. _config-env:

Scope `env`
-----------

The ``env`` scope allows the definition one or more variable that will be exported in the environment where the
workflow tasks will be executed.

Simply prefix your variable names with the ``env`` scope or surround them by curly brackets, as shown below::

   env.ALPHA = 'some value'
   env.BETA = "$HOME/some/path"

   env {
        DELTA = 'one more'
        GAMMA = "/my/path:$PATH"
   }

.. note::
  In the above example, variables like ``$HOME`` and ``$PATH`` are evaluated when the workflow is launched. If
  you want these variables to be evaluated during task execution, escape them with ``\$``. This difference is important
  for variables like ``$PATH``, which may be different in the workflow environment versus the task environment.

.. warning::
  The ``env`` scope provides environment variables to *tasks*, not Nextflow itself. Nextflow environment variables
  such as ``NXF_VER`` should be set in the environment in which Nextflow is launched.


.. _config-executor:

Scope `executor`
----------------

The ``executor`` configuration scope allows you to set the optional executor settings, listed in the following table.

===================== =====================
Name                  Description
===================== =====================
name                  The name of the executor to be used (default: ``local``).
queueSize             The number of tasks the executor will handle in a parallel manner (default: ``100``).
submitRateLimit       Determines the max rate of job submission per time unit, for example ``'10sec'`` (10 jobs per second) or ``'50/2min'`` (50 jobs every 2 minutes) (default: unlimited).
pollInterval          Determines how often to check for process termination. Default varies for each executor.
dumpInterval          Determines how often to log the executor status (default: ``5min``).
queueStatInterval     Determines how often to fetch the queue status from the scheduler (default: ``1min``). Used only by grid executors.
exitReadTimeout       Determines how long to wait before returning an error status when a process is terminated but the ``.exitcode`` file does not exist or is empty (default: ``270 sec``). Used only by grid executors.
killBatchSize         Determines the number of jobs that can be killed in a single command execution (default: ``100``).
perJobMemLimit        Specifies Platform LSF *per-job* memory limit mode. See :ref:`lsf-executor`.
perTaskReserve        Specifies Platform LSF *per-task* memory reserve mode. See :ref:`lsf-executor`.
jobName               Determines the name of jobs submitted to the underlying cluster executor e.g. ``executor.jobName = { "$task.name - $task.hash" }``. Make sure the resulting job name matches the validation constraints of the underlying batch scheduler.
cpus                  The maximum number of CPUs made available by the underlying system. Used only by the ``local`` executor.
memory                The maximum amount of memory made available by the underlying system. Used only by the ``local`` executor.
retry.delay           Delay when retrying failed job submissions (default: ``500ms``). NOTE: used only by grid executors (requires ``22.03.0-edge`` or later).
retry.maxDelay        Max delay when retrying failed job submissions (default: ``30s``). NOTE: used only by grid executors (requires ``22.03.0-edge`` or later).
retry.jitter          Jitter value when retrying failed job submissions (default: ``0.25``). NOTE: used only by grid executors (requires ``22.03.0-edge`` or later).
retry.maxAttempts     Max attempts when retrying failed job submissions (default: ``3``). NOTE: used only by grid executors (requires ``22.03.0-edge`` or later).
retry.reason          Regex pattern that when verified cause a failed submit operation to be re-tried (default: ``Socket timed out``). NOTE: used only by grid executors (requires ``22.03.0-edge`` or later).
===================== =====================

The executor settings can be defined as shown below::

    executor {
        name = 'sge'
        queueSize = 200
        pollInterval = '30 sec'
    }

When using two (or more) different executors in your pipeline, you can specify their settings separately by prefixing
the executor name with the symbol ``$`` and using it as special scope identifier. For example::

  executor {
    $sge {
        queueSize = 100
        pollInterval = '30sec'
    }

    $local {
        cpus = 8
        memory = '32 GB'
    }
  }

The above configuration example can be rewritten using the dot notation as shown below::

    executor.$sge.queueSize = 100
    executor.$sge.pollInterval = '30sec'
    executor.$local.cpus = 8
    executor.$local.memory = '32 GB'


.. _config-k8s:

Scope `k8s`
-----------

The ``k8s`` scope allows the definition of the configuration settings that control the deployment and execution of
workflow applications in a Kubernetes cluster.

The following settings are available:

=================== ================
Name                Description
=================== ================
autoMountHostPaths  Automatically mounts host paths in the job pods. Only for development purpose when using a single node cluster (default: ``false``).
context             Defines the Kubernetes `configuration context name <https://kubernetes.io/docs/tasks/access-application-cluster/configure-access-multiple-clusters/>`_ to use.
namespace           Defines the Kubernetes namespace to use (default: ``default``).
serviceAccount      Defines the Kubernetes `service account name <https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/>`_ to use.
launchDir           Defines the path where the workflow is launched and the user data is stored. This must be a path in a shared K8s persistent volume (default: ``<volume-claim-mount-path>/<user-name>``.
workDir             Defines the path where the workflow temporary data is stored. This must be a path in a shared K8s persistent volume (default:``<user-dir>/work``).
projectDir          Defines the path where Nextflow projects are downloaded. This must be a path in a shared K8s persistent volume (default: ``<volume-claim-mount-path>/projects``).
pod                 Allows the definition of one or more pod configuration options such as environment variables, config maps, secrets, etc. It allows the same settings as the :ref:`process-pod` process directive.
pullPolicy          Defines the strategy to be used to pull the container image e.g. ``pullPolicy: 'Always'``.
runAsUser           Defines the user ID to be used to run the containers. Shortcut for the ``securityContext`` option.
securityContext     Defines the `security context <https://kubernetes.io/docs/tasks/configure-pod-container/security-context/>`_ for all pods.
storageClaimName    The name of the persistent volume claim where store workflow result data.
storageMountPath    The path location used to mount the persistent volume claim (default: ``/workspace``).
storageSubPath      The path in the persistent volume to be mounted (default: root).
computeResourceType Define whether use Kubernetes ``Pod`` or ``Job`` resource type to carry out Nextflow tasks (default: ``Pod``).
fetchNodeName       If you trace the hostname, activate this option (default: ``false``, requires version ``22.05.0-edge`` or later).
volumeClaims        (deprecated)
maxErrorRetry       Defines the Kubernetes API max request retries (default is set to 4)
httpReadTimeout     Defines the Kubernetes client request HTTP connection read timeout e.g. ``'60s'`` (requires version ``22.10.0`` or later).
httpConnectTimeout  Defines the Kubernetes client request HTTP connection timeout e.g. ``'60s'`` (requires version ``22.10.0`` or later).
=================== ================

See the :ref:`k8s-page` documentation for more details.


.. _config-mail:

Scope `mail`
------------

The ``mail`` scope allows you to define the mail server configuration settings needed to send email messages.

================== ================
Name                Description
================== ================
from                Default email sender address.
smtp.host           Host name of the mail server.
smtp.port           Port number of the mail server.
smtp.user           User name to connect to  the mail server.
smtp.password       User password to connect to the mail server.
smtp.proxy.host     Host name of an HTTP web proxy server that will be used for connections to the mail server.
smtp.proxy.port     Port number for the HTTP web proxy server.
smtp.*              Any SMTP configuration property supported by the Java Mail API (see link below).
debug               When ``true`` enables Java Mail logging for debugging purpose.
================== ================

.. note:: Nextflow relies on the `Java Mail API <https://javaee.github.io/javamail/>`_ to send email messages.
  Advanced mail configuration can be provided by using any SMTP configuration property supported by the Java Mail API.
  See the `table of available properties at this link <https://javaee.github.io/javamail/docs/api/com/sun/mail/smtp/package-summary.html#properties>`_.

For example, the following snippet shows how to configure Nextflow to send emails through the
`AWS Simple Email Service <https://aws.amazon.com/ses/>`_::

    mail {
        smtp.host = 'email-smtp.us-east-1.amazonaws.com'
        smtp.port = 587
        smtp.user = '<Your AWS SES access key>'
        smtp.password = '<Your AWS SES secret key>'
        smtp.auth = true
        smtp.starttls.enable = true
        smtp.starttls.required = true
    }


.. _config-manifest:

Scope `manifest`
----------------

The ``manifest`` configuration scope allows you to define some meta-data information needed when publishing your pipeline project on GitHub, BitBucket or GitLab, or when running your pipeline.

The following settings are available:

================== ================
Name                Description
================== ================
author              Project author name (use a comma to separate multiple names).
defaultBranch       Git repository default branch (default: ``master``).
recurseSubmodules   Turn this flag to ``true`` to pull submodules recursively from the Git repository
description         Free text describing the workflow project.
doi                 Project related publication DOI identifier.
homePage            Project home page URL.
mainScript          Project main script (default: ``main.nf``).
name                Project short name.
nextflowVersion     Minimum required Nextflow version.
version             Project version number.
================== ================

The above options can be used by prefixing them with the ``manifest`` scope or surrounding them by curly
brackets. For example::

    manifest {
        homePage = 'http://foo.com'
        description = 'Pipeline does this and that'
        mainScript = 'foo.nf'
        version = '1.0.0'
    }

To learn how to publish your pipeline on GitHub, BitBucket or GitLab code repositories read :ref:`sharing-page`
documentation page.

Nextflow version
^^^^^^^^^^^^^^^^

The ``nextflowVersion`` setting allows you to specify a minimum required version to run the pipeline.
This may be useful to ensure that a specific version is used::

    nextflowVersion = '1.2.3'        // exact match
    nextflowVersion = '1.2+'         // 1.2 or later (excluding 2 and later)
    nextflowVersion = '>=1.2'        // 1.2 or later
    nextflowVersion = '>=1.2, <=1.5' // any version in the 1.2 .. 1.5 range
    nextflowVersion = '!>=1.2'       // with ! prefix, stop execution if current version
                                        does not match required version.


.. _config-notification:

Scope `notification`
--------------------

The ``notification`` scope allows you to define the automatic sending of a notification email message
when the workflow execution terminates.

================== ================
Name                Description
================== ================
enabled             Enables the sending of a notification message when the workflow execution completes.
to                  Recipient address for the notification email. Multiple addresses can be specified separating them with a comma.
from                Sender address for the notification email message.
template            Path of a template file which provides the content of the notification message.
binding             An associative array modelling the variables in the template file.
================== ================

The notification message is sent my using the STMP server defined in the configuration :ref:`mail scope<config-mail>`.

If no mail configuration is provided, it tries to send the notification message by using the external mail command
eventually provided by the underlying system (eg. ``sendmail`` or ``mail``).


.. _config-params:

Scope `params`
--------------

The ``params`` scope allows you to define parameters that will be accessible in the pipeline script. Simply prefix the
parameter names with the ``params`` scope or surround them by curly brackets, as shown below::

    params.custom_param = 123
    params.another_param = 'string value .. '

    params {
        alpha_1 = true
        beta_2 = 'another string ..'
    }


.. _config-podman:

Scope `podman`
--------------

The ``podman`` configuration scope controls how `Podman <https://podman.io/>`_ containers are executed by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Podman execution (default: ``false``).
envWhitelist        Comma separated list of environment variable names to be included in the container environment.
temp                Mounts a path of your choice as the ``/tmp`` directory in the container. Use the special value ``auto`` to create a temporary directory each time a container is created.
remove              Clean-up the container after the execution (default: ``true``).
runOptions          This attribute can be used to provide any extra command line options supported by the ``podman run`` command.
registry            The registry from where container images are pulled. It should be only used to specify a private registry server. It should NOT include the protocol prefix i.e. ``http://``.
engineOptions       This attribute can be used to provide any option supported by the Docker engine i.e. ``podman [OPTIONS]``.
mountFlags          Add the specified flags to the volume mounts e.g. `mountFlags = 'ro,Z'`
================== ================

The above options can be used by prefixing them with the ``podman`` scope or surrounding them by curly
brackets, as shown below::

    process.container = 'nextflow/examples'

    podman {
        enabled = true
        temp = 'auto'
    }

Read :ref:`container-podman` page to learn more about how to use Podman containers with Nextflow.


.. _config-process:

Scope `process`
---------------

The ``process`` configuration scope allows you to provide the default configuration for the processes in your pipeline.

You can specify here any property described in the :ref:`process directive<process-directives>` and the executor sections.
For examples::

    process {
        executor = 'sge'
        queue = 'long'
        clusterOptions = '-pe smp 10 -l virtual_free=64G,h_rt=30:00:00'
    }

By using this configuration all processes in your pipeline will be executed through the SGE cluster, with the specified
settings.


.. _config-process-selectors:

Process selectors
^^^^^^^^^^^^^^^^^

The ``withLabel`` selectors allow the configuration of all processes annotated with a :ref:`process-label` directive as
shown below::

    process {
        withLabel: big_mem {
            cpus = 16
            memory = 64.GB
            queue = 'long'
        }
    }

The above configuration example assigns 16 cpus, 64 Gb of memory and the ``long`` queue to all processes annotated
with the ``big_mem`` label.

In the same manner, the ``withName`` selector allows the configuration of a specific process in your pipeline by its name.
For example::

    process {
        withName: hello {
            cpus = 4
            memory = 8.GB
            queue = 'short'
        }
    }

.. tip::
  Label and process names do not need to be enclosed with quotes, provided the name
  does not include special characters (``-``, ``!``, etc) and is not a keyword or a built-in type identifier.
  When in doubt, you can enclose the label name or process name with single or double quotes.


.. _config-selector-expressions:

Selector expressions
^^^^^^^^^^^^^^^^^^^^

Both label and process name selectors allow the use of a regular expression in order to apply the same configuration
to all processes matching the specified pattern condition. For example::

    process {
        withLabel: 'foo|bar' {
            cpus = 2
            memory = 4.GB
        }
    }

The above configuration snippet sets 2 cpus and 4 GB of memory to the processes annotated with with a label ``foo``
and ``bar``.

A process selector can be negated prefixing it with the special character ``!``. For example::

    process {
        withLabel: 'foo' { cpus = 2 }
        withLabel: '!foo' { cpus = 4 }
        withName: '!align.*' { queue = 'long' }
    }

The above configuration snippet sets 2 cpus for the processes annotated with the ``foo`` label and 4 cpus to all processes
*not* annotated with that label. Finally it sets the use of ``long`` queue to all process whose name does *not* start
with ``align``.


.. _config-selector-priority:

Selector priority
^^^^^^^^^^^^^^^^^

When mixing generic process configuration and selectors the following priority rules are applied (from lower to higher):

1. Process generic configuration.
2. Process specific directive defined in the workflow script.
3. ``withLabel`` selector definition.
4. ``withName`` selector definition.

For example::

    process {
        cpus = 4
        withLabel: foo { cpus = 8 }
        withName: bar { cpus = 32 }
    }

Using the above configuration snippet, all workflow processes use 4 cpus if not otherwise specified in the workflow
script. Moreover processes annotated with the ``foo`` label use 8 cpus. Finally the process named ``bar``
uses 32 cpus.


.. _config-report:

Scope `report`
--------------

The ``report`` scope allows you to define configuration setting of the workflow :ref:`execution-report`.

================== ================
Name                Description
================== ================
enabled             If ``true`` it create the workflow execution report.
file                The path of the created execution report file (default: ``report.html``).
overwrite           When ``true`` overwrites existing report file instead of rolling it.
================== ================


.. _config-shifter:

Scope `shifter`
-------------------

The ``shifter`` configuration scope controls how `Shifter <https://docs.nersc.gov/programming/shifter/overview/>`_ containers are executed
by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Shifter execution (default: ``false``).
================== ================

Read :ref:`container-shifter` page to learn more about how to use Shifter containers with Nextflow.


.. _config-singularity:

Scope `singularity`
-------------------

The ``singularity`` configuration scope controls how `Singularity <https://sylabs.io/singularity/>`_ containers are executed
by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Singularity execution (default: ``false``).
engineOptions       This attribute can be used to provide any option supported by the Singularity engine i.e. ``singularity [OPTIONS]``.
envWhitelist        Comma separated list of environment variable names to be included in the container environment.
runOptions          This attribute can be used to provide any extra command line options supported by the ``singularity exec``.
noHttps             Turn this flag to ``true`` to pull the Singularity image with http protocol (default: ``false``).
autoMounts          When ``true`` Nextflow automatically mounts host paths in the executed container. It requires the `user bind control` feature enabled in your Singularity installation (default: ``false``).
cacheDir            The directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible to all compute nodes.
pullTimeout         The amount of time the Singularity pull can last, exceeding which the process is terminated (default: ``20 min``).
================== ================

Read :ref:`container-singularity` page to learn more about how to use Singularity containers with Nextflow.


.. _config-timeline:

Scope `timeline`
----------------

The ``timeline`` scope allows you to enable/disable the processes execution timeline report generated by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             When ``true`` turns on the generation of the timeline report file (default: ``false``).
file                Timeline file name (default: ``timeline.html``).
overwrite           When ``true`` overwrites an existing timeline file instead of rolling it.
================== ================


.. _config-tower:

Scope `tower`
-------------

The ``tower`` configuration scope controls the settings for the `Nextflow Tower <https://tower.nf>`_ monitoring and tracing service.

The following settings are available:

================== ================
Name                Description
================== ================
enabled            When ``true`` Nextflow sends the workflow tracing and execution metrics to the Nextflow Tower service (default: ``false``).
accessToken        The unique access token specific to your account on an instance of Tower.
endpoint           The endpoint of your Tower deployment (default: ``https://tower.nf``).
workspaceId        The ID of the Tower workspace where the run should be added (default: the launching user personal workspace).
================== ================

The above options can be used by prefixing them with the ``tower`` scope or surrounding them by curly
brackets, as shown below::

    tower {
      enabled = true
      accessToken = '<YOUR TOKEN>'
      workspaceId = '<YOUR WORKSPACE ID>'
    }

.. tip::
  Your ``accessToken`` can be obtained from your Tower instance in the `Tokens page <https://tower.nf/tokens>`.

.. tip::
  The Tower workspace ID can also be specified using the environment variable ``TOWER_WORKSPACE_ID`` (config file has priority over the environment variable).


.. _config-trace:

Scope `trace`
-------------

The ``trace`` scope allows you to control the layout of the execution trace file generated by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             When ``true`` turns on the generation of the execution trace report file (default: ``false``).
fields              Comma separated list of fields to be included in the report. The available fields are listed at :ref:`this page <trace-fields>`
file                Trace file name (default: ``trace.txt``).
sep                 Character used to separate values in each row (default: ``\t``).
raw                 When ``true`` turns on raw number report generation i.e. date and time are reported as milliseconds and memory as number of bytes
overwrite           When ``true`` overwrites an existing trace file instead of rolling it.
================== ================

The above options can be used by prefixing them with the ``trace`` scope or surrounding them by curly
brackets. For example::

    trace {
        enabled = true
        file = 'pipeline_trace.txt'
        fields = 'task_id,name,status,exit,realtime,%cpu,rss'
    }

To learn more about the execution report that can be generated by Nextflow read :ref:`trace-report` documentation page.


.. _config-weblog:

Scope `weblog`
--------------

The ``weblog`` scope allows you to send detailed :ref:`trace scope<trace-fields>` information as HTTP POST request to a webserver, shipped as a JSON object.

Detailed information about the JSON fields can be found in the :ref:`weblog description<weblog-service>`.

================== ================
Name                Description
================== ================
enabled             If ``true`` it will send HTTP POST requests to a given url.
url                The url where to send HTTP POST requests (default: ``http:localhost``).
================== ================


.. _config-miscellaneous:

Miscellaneous
-------------

There are additional variables that can be defined within a configuration file that do not have a dedicated scope.

These are defined alongside other scopes, but the option is assigned as typically variable.

================== ================
Name                Description
================== ================
cleanup             If ``true``, on a successful completion of a run all files in *work* directory are automatically deleted.
================== ================

.. warning::
    The use of the ``cleanup`` option will prevent the use of the *resume* feature on subsequent executions of that pipeline run.
    Also, be aware that deleting all scratch files can take a lot of time, especially when using a shared file system or remote cloud storage.


.. _config-profiles:

Config profiles
===============

Configuration files can contain the definition of one or more *profiles*. A profile is a set of configuration attributes
that can be activated/chosen when launching a pipeline execution by using the ``-profile`` command line option.

Configuration profiles are defined by using the special scope ``profiles`` which group the attributes that belong
to the same profile using a common prefix. For example::

    profiles {

        standard {
            process.executor = 'local'
        }

        cluster {
            process.executor = 'sge'
            process.queue = 'long'
            process.memory = '10GB'
        }

        cloud {
            process.executor = 'cirrus'
            process.container = 'cbcrg/imagex'
            docker.enabled = true
        }

    }

This configuration defines three different profiles: ``standard``, ``cluster`` and ``cloud`` that set different process
configuration strategies depending on the target runtime platform. By convention the ``standard`` profile is implicitly used
when no other profile is specified by the user.

.. tip::
    Multiple configuration profiles can be specified by separating the profile names
    with a comma, for example::

        nextflow run <your script> -profile standard,cloud

.. danger::
    When using the ``profiles`` feature in your config file, do NOT set attributes in the same scope both
    inside and outside a ``profiles`` context. For example::

        process.cpus = 1

        profiles {
          foo {
            process.memory = '2 GB'
          }

          bar {
            process.memory = '4 GB'
          }
        }

    In the above example, the ``process.cpus`` attribute is not correctly applied because the ``process`` scope is also
    used in the ``foo`` and ``bar`` profiles.


.. _config-env-vars:

Environment variables
=====================

The following environment variables control the configuration of the Nextflow runtime and
the underlying Java virtual machine.

=============================== ================
Name                            Description
=============================== ================
NXF_ANSI_LOG                    Enables/disables ANSI console output (default ``true`` when ANSI terminal is detected).
NXF_ANSI_SUMMARY                Enables/disables ANSI completion summary: `true|false` (default: print summary if execution last more than 1 minute).
NXF_ASSETS                      Defines the directory where downloaded pipeline repositories are stored (default: ``$NXF_HOME/assets``)
NXF_CHARLIECLOUD_CACHEDIR       Directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.
NXF_CLASSPATH                   Allows the extension of the Java runtime classpath with extra JAR files or class folders.
NXF_CLOUD_DRIVER                Defines the default cloud driver to be used if not specified in the config file or as command line option, either ``aws`` or ``google``.
NXF_CONDA_CACHEDIR              Directory where Conda environments are store. When using a computing cluster it must be a shared folder accessible from all compute nodes.
NXF_CONDA_ENABLED               Enable the use of Conda recipes defined by using the :ref:process-conda directive. (default: ``false``, requires version ``22.08.0-edge`` or later).
NXF_DEBUG                       Defines scripts debugging level: ``1`` dump task environment variables in the task log file; ``2`` enables command script execution tracing; ``3`` enables command wrapper execution tracing.
NXF_DEFAULT_DSL                 Defines the DSL version version that should be used in not specified otherwise in the script of config file (default: ``2``, requires version ``22.03.0-edge`` or later)
NXF_DISABLE_JOBS_CANCELLATION   Disables the cancellation of child jobs on workflow execution termination (requires version ``21.12.0-edge`` or later).
NXF_ENABLE_STRICT               Enable Nextflow *strict* execution mode (default: ``false``, requires version ``22.05.0-edge`` or later)
NXF_ENABLE_SECRETS              Enable Nextflow secrets features (default: ``true``, requires version ``21.09.0-edge`` or later)
NXF_EXECUTOR                    Defines the default process executor e.g. `sge`
NXF_GRAB                        Provides extra runtime dependencies downloaded from a Maven repository service [DEPRECATED]
NXF_HOME                        Nextflow home directory (default: ``$HOME/.nextflow``).
NXF_JAVA_HOME                   Defines the path location of the Java VM installation used to run Nextflow. This variable overrides the ``JAVA_HOME`` variable if defined.
NXF_JVM_ARGS                    Allows the setting Java VM options. This is similar to ``NXF_OPTS`` however it's only applied the JVM running Nextflow and not to any java pre-launching commands (requires ``21.12.1-edge`` or later).
NXF_OFFLINE                     When ``true`` disables the project automatic download and update from remote repositories (default: ``false``).
NXF_OPTS                        Provides extra options for the Java and Nextflow runtime. It must be a blank separated list of ``-Dkey[=value]`` properties.
NXF_ORG                         Default `organization` prefix when looking for a hosted repository (default: ``nextflow-io``).
NXF_PARAMS_FILE                 Defines the path location of the pipeline parameters file (requires version ``20.10.0`` or later).
NXF_PID_FILE                    Name of the file where the process PID is saved when Nextflow is launched in background.
NXF_SCM_FILE                    Defines the path location of the SCM config file (requires version ``20.10.0`` or later).
NXF_SINGULARITY_CACHEDIR        Directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.
NXF_SINGULARITY_LIBRARYDIR      Directory where remote Singularity images are retrieved. It should be a directory accessible to all compute nodes (requires: ``21.09.0-edge`` or later).
NXF_TEMP                        Directory where temporary files are stored
NXF_VER                         Defines what version of Nextflow to use.
NXF_WORK                        Directory where working files are stored (usually your *scratch* directory)
JAVA_HOME                       Defines the path location of the Java VM installation used to run Nextflow.
JAVA_CMD                        Defines the path location of the Java binary command used to launch Nextflow.
HTTP_PROXY                      Defines the HTTP proxy server. As of version ``21.06.0-edge``, proxy authentication is supported providing the credentials in the proxy URL e.g. ``http://user:password@proxy-host.com:port``.
HTTPS_PROXY                     Defines the HTTPS proxy server. As of version ``21.06.0-edge``, proxy authentication is supported providing the credentials in the proxy URL e.g. ``https://user:password@proxy-host.com:port``.
FTP_PROXY                       Defines the FTP proxy server. Proxy authentication is supported providing the credentials in the proxy URL e.g. ``ftp://user:password@proxy-host.com:port``. FTP proxy support requires version ``21.06.0-edge`` or later.
NO_PROXY                        Defines one or more host names that should not use the proxy server. Separate multiple names using a comma character.
=============================== ================
