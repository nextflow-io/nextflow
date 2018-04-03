.. _config-page:

*************
Configuration
*************

Configuration file
==================

When a pipeline script is launched Nextflow looks for a file named ``nextflow.config`` in the current directory and
in the script base directory (if it is not the same as the current directory). Finally it checks for the file
``$HOME/.nextflow/config``.

When more than one on the above files exist they are merged, so that the settings in the first override the same ones
that may appear in the second one, and so on.

The default config file search mechanism can be extended proving an extra configuration file by using the command line
option ``-c <config file>``.

.. note:: It's worth noting that by doing this, the files ``nextflow.config`` and ``$HOME/.nextflow/config`` are not
  ignored and they are merged as explained above.

.. tip:: If you want to ignore any default configuration files and use only the custom one use the command line option
  ``-C <config file>``.

Config syntax
--------------

A Nextflow configuration file is a simple text file containing a set of properties defined using the syntax::

  name = value

Please note, string values need to be wrapped in quotation characters while numbers and boolean values (``true``, ``false``) do not.
Also note that values are typed, meaning for example that, ``1`` is different from ``'1'``, since the first is interpreted
as the number one, while the latter is interpreted as a string value.


Config Variables
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

.. warning:: If you are accessing an environment variable that may not exist in the system, your property may contain
    an undefined value. You can avoid this by using a conditional expression in your property definition as shown below.

::

    mySafeProperty = "${MY_FANCY_VARIABLE?:''}"


Config comments
------------------

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



Scope `env`
-----------

The ``env`` scope allows you to define one or more environment variables that will be exported to the system environment
where pipeline processes need to be executed.

Simply prefix your variable names with the ``env`` scope or surround them by curly brackets, as shown below::

   env.ALPHA = 'some value'
   env.BETA = "$HOME/some/path"

   env {
        DELTA = 'one more'
        GAMMA = "/my/path:$PATH"
   }


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



.. _config-process:

Scope `process`
---------------

The ``process`` configuration scope allows you to provide the default configuration for the processes in your pipeline.

You can specify here any property described in the :ref:`process directive<process-directives>` and the executor sections.
For examples::

  process {
    executor='sge'
    queue='long'
    clusterOptions = '-pe smp 10 -l virtual_free=64G,h_rt=30:00:00'
  }


By using this configuration all processes in your pipeline will be executed through the SGE cluster, with the specified
settings.

.. _config-process-selectors:

Process selectors
^^^^^^^^^^^^^^^^^

The ``withLabel`` selectors allow the configuration of all processes annotated with a :ref:`process-label` directive
shown below::

    process {
        withLabel: big_mem {
            cpus = 16
            memory = 64.GB
            queue = 'long'
        }
    }

The above configuration example assigns 16 cpus, 64 Gb of memory and the ``long`` queue to all processes annotated
with a ``big_mem`` label.


In the same manner, the ``withName`` selector allows the configuration of a specific process in your pipeline by its name.
For example::

    process {
        withName: hello {
            cpus = 4
            memory = 8.GB
            queue = 'short'
        }
    }

.. tip:: Either label and process names do not need to be enclosed with quote characters, provided the name
  does include special characters (e.g. ``-``, ``!``, etc) or it's not a keyword or a built-in type identifier.
  In case of doubt, you can enclose the label name or the process name with single or double quote characters.

.. _config-selector-expressions:

Selector expressions
^^^^^^^^^^^^^^^^^^^^

Both label and process name selectors allow the use of a regular expression to apply the same configuration
to all processes matching the selection pattern condition. For example::

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

The above configuration snippet sets 2 cpus to processes annotated with the ``foo`` label and 4 cpus to all processes
*not* annotated with that label. Finally it sets the use of ``long`` queue to all process whose name does *not* start
with ``align``.

.. _config-selectors-priority:

Selectors priority
^^^^^^^^^^^^^^^^^^

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

Using the above configuration snippet, all workflow processes use 4 cpus in not otherwise specified in the workflow
script. Moreover processes annotated with the ``foo`` label use 8 cpus. Finally the process named ``bar``
uses 32 cpus.


.. _config-executor:

Scope `executor`
----------------

The ``executor`` configuration scope allows you to set the optional executor settings, listed in the following table.

===================== =====================
Name                  Description
===================== =====================
name                  The name of the executor to be used e.g. ``local``, ``sge``, etc.
queueSize             The number of tasks the executor will handle in a parallel manner (default: ``100``).
pollInterval          Determines how often a poll occurs to check for a process termination.
dumpInterval          Determines how often the executor status is written in the application log file (default: ``5min``).
queueStatInterval     Determines how often the queue status is fetched from the cluster system. This setting is used only by grid executors (default: ``1min``).
exitReadTimeout       Determines how long the executor waits before return an error status when a process is terminated but the `exit` file does not exist or it is empty. This setting is used only by grid executors (default: ``270 sec``).
killBatchSize         Determines the number of jobs that can be `killed` in a single command execution (default: ``100``).
perJobMemLimit        Specifies Platform LSF *per-job* memory limit mode. See :ref:`lsf-executor`.
jobName               Determines the name of jobs submitted to the underlying cluster executor e.g. ``executor.jobName = { "$task.name - $task.hash" }`` .
cpus                  The maximum number of CPUs made available by the underlying system (only used by the ``local`` executor).
memory                The maximum amount of memory made available by the underlying system (only used by the ``local`` executor).
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

.. _config-docker:

Scope `docker`
--------------

The ``docker`` configuration scope controls how `Docker <http://www.docker.io>`_ containers are executed by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Docker execution (default: ``false``).
legacy              Uses command line options removed since version 1.10.x (default: ``false``).
sudo                Executes Docker run command as ``sudo`` (default: ``false``).
tty                 Allocates a pseudo-tty (default: ``false``).
temp                Mounts a path of your choice as the ``/tmp`` directory in the container. Use the special value ``auto`` to create a temporary directory each time a container is created.
remove              Clean-up the container after the execution (default: ``true``). For details see: http://docs.docker.com/reference/run/#clean-up-rm .
runOptions          This attribute can be used to provide any extra command line options supported by the ``docker run`` command. For details see: http://docs.docker.com/reference/run .
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



Read :ref:`docker-page` page to lean more how use Docker containers with Nextflow.


.. _config-singularity:

Scope `singularity`
-------------------

The ``singularity`` configuration scope controls how `Singularity <http://singularity.lbl.gov>`_ containers are executed
by Nextflow.

The following settings are available:

================== ================
Name                Description
================== ================
enabled             Turn this flag to ``true`` to enable Singularity execution (default: ``false``).
engineOptions       This attribute can be used to provide any option supported by the Singularity engine i.e. ``singularity [OPTIONS]``.
runOptions          This attribute can be used to provide any extra command line options supported by the ``singularity exec``.
autoMounts          When ``true`` Nextflow automatically mounts host paths in the executed contained. It requires the `user bind control` feature enabled in your Singularity installation (default: ``false``).
cacheDir            The directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible to all computing nodes.
================== ================


Read :ref:`singularity-page` page to lean more how use Singularity containers with Nextflow.

.. _config-manifest:

Scope `manifest`
----------------

The ``manifest`` configuration scope allows you to define some meta-data information needed when publishing your
pipeline project on GitHub, BitBucket or GitLab.

The following settings are available:

================== ================
Name                Description
================== ================
author              Project author name (use a comma to separate multiple names).
homePage            Project home page URL
description         Free text describing the pipeline project
mainScript          Pipeline main script (default: ``main.nf``)
defaultBranch       Git repository default branch (default: ``master``)
================== ================

The above options can be used by prefixing them with the ``manifest`` scope or surrounding them by curly
brackets. For example::

    manifest {
        homePage = 'http://foo.com'
        description = 'Pipeline does this and that'
        mainScript = 'foo.nf'
    }


To learn how to publish your pipeline on GitHub, BitBucket or GitLab code repositories read :ref:`sharing-page`
documentation page.

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
================== ================

The above options can be used by prefixing them with the ``trace`` scope or surrounding them by curly
brackets. For example::

    trace {
        enabled = true
        file = 'pipeline_trace.txt'
        fields = 'task_id,name,status,exit,realtime,%cpu,rss'
    }


To learn more about the execution report that can be generated by Nextflow read :ref:`trace-report` documentation page.

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

Click the following link to lean more about `AWS Security Credentials <http://docs.aws.amazon.com/general/latest/gr/aws-security-credentials.html>`_.

Advanced client configuration options can be set by using the ``client`` attribute. The following properties can be used:

=========================== ================
Name                        Description
=========================== ================
connectionTimeout           The amount of time to wait (in milliseconds) when initially establishing a connection before giving up and timing out.
endpoint                    The AWS S3 API entry point e.g. `s3-us-west-1.amazonaws.com`.
maxConnections              The maximum number of allowed open HTTP connections.
maxErrorRetry               The maximum number of retry attempts for failed retryable requests.
protocol                    The protocol (i.e. HTTP or HTTPS) to use when connecting to AWS.
proxyHost                   The proxy host to connect through.
proxyPort                   The port on the proxy host to connect through.
proxyUsername               The user name to use when connecting through a proxy.
proxyPassword               The password to use when connecting through a proxy.
signerOverride              The name of the signature algorithm to use for signing requests made by the client.
socketSendBufferSizeHint    The Size hint (in bytes) for the low level TCP send buffer.
socketRecvBufferSizeHint    The Size hint (in bytes) for the low level TCP receive buffer.
socketTimeout               The amount of time to wait (in milliseconds) for data to be transferred over an established, open connection before the connection is timed out.
storageEncryption           The S3 server side encryption to be used when saving objects on S3 (currently only AES256 is supported)
userAgent                   The HTTP user agent header passed with all HTTP requests.
uploadMaxThreads            The maximum number of threads used for multipart upload.
uploadChunkSize             The size of a single part in a multipart upload (default: `10 MB`).
uploadStorageClass          The S3 storage class applied to stored objects, either `STANDARD` or `REDUCED_REDUNDANCY` (default: `STANDARD`).
uploadMaxAttempts           The maximum number of upload attempts after which a multipart upload returns an error (default: `5`).
uploadRetrySleep            The time to wait after a failed upload attempt to retry the part upload (default: `100ms`).
=========================== ================

For example::

    aws {
        client {
            maxConnections = 20
            connectionTimeout = 10000
            uploadStorageClass = 'REDUCED_REDUNDANCY'
            storageEncryption = 'AES256'
        }
    }


.. _config-cloud:

Scope `cloud`
-------------

The ``cloud`` scope allows you to define the settings of the computing cluster that can be deployed in the cloud
by Nextflow.

The following settings are available:

=========================== ================
Name                        Description
=========================== ================
bootStorageSize             Boot storage volume size e.g. ``10 GB``.
imageId                     Identifier of the virtual machine(s) to launch e.g. ``ami-43f49030``.
instanceRole                IAM role granting required permissions and authorizations in the launched instances.
                            When specifying an IAM role no access/security keys are installed in the cluster deployed in the cloud.
instanceType                Type of the virtual machine(s) to launch e.g. ``m4.xlarge``.
instanceStorageMount        Ephemeral instance storage mount path e.g. ``/mnt/scratch``.
instanceStorageDevice       Ephemeral instance storage device name e.g. ``/dev/xvdc`` (optional).
keyName                     SSH access key name given by the cloud provider.
keyHash                     SSH access public key hash string.
keyFile                     SSH access public key file path.
securityGroup               Identifier of the security group to be applied e.g. ``sg-df72b9ba``.
sharedStorageId             Identifier of the shared file system instance e.g. ``fs-1803efd1``.
sharedStorageMount          Mount path of the shared file system e.g. ``/mnt/efs``.
subnetId                    Identifier of the VPC subnet to be applied e.g. ``subnet-05222a43``.
spotPrice                   Price bid for spot/preemptive instances.
userName                    SSH access user name (don't specify it to use the image default user name).
autoscale                   See below.
=========================== ================

The autoscale configuration group provides the following settings:

=========================== ================
Name                        Description
=========================== ================
enabled                     Enable cluster auto-scaling.
terminateWhenIdle           Enable cluster automatic scale-down i.e. instance terminations when idle (default: ``false``).
idleTimeout                 Amount of time in idle state after which an instance is candidate to be terminated (default: ``5 min``).
starvingTimeout             Amount of time after which one ore more tasks pending for execution trigger an auto-scale request (default: ``5 min``).
minInstances                Minimum number of instances in the cluster.
maxInstances                Maximum number of instances in the cluster.
imageId                     Identifier of the virtual machine(s) to launch when new instances are added to the cluster.
instanceType                Type of the virtual machine(s) to launch when new instances are added to the cluster.
spotPrice                   Price bid for spot/preemptive instances launched while auto-scaling the cluster.
=========================== ================

.. _config-k8s:

Scope `k8s`
-----------

The ``k8s`` scope allows the definition of the configuration settings that control the deployment and execution of
workflow applications in a Kubernetes cluster.

The following settings are available:

================== ================
Name                Description
================== ================
context             Defines the Kubernetes `configuration context name <https://kubernetes.io/docs/tasks/access-application-cluster/configure-access-multiple-clusters/>`_ to use.
namespace           Defines the Kubernetes namespace to use (default: ``default``).
serviceAccount      Defines the Kubernetes `service account name <https://kubernetes.io/docs/tasks/configure-pod-container/configure-service-account/>`_ to use.
userDir             Defines the path where the workflow is launched and the user data is stored. This must be a path in a shared K8s persistent volume (default: ``<volume-claim-mount-path>/<user-name>``.
workDir             Defines the path where the workflow temporary data is stored. This must be a path in a shared K8s persistent volume (default:``<user-dir>/work``).
projectDir          Defines the path where Nextflow projects are downloaded. This must be a path in a shared K8s persistent volume (default: ``<volume-claim-mount-path>/projects``).
volumeClaims        Configures one or more persistent volume claims in the execution environment. See below for details.
================== ================

Volume claims need to be defined as a named object specifying the ``mountPath`` as a nested object::

    k8s {
      volumeClaims = [ 'your-pvc-name': [mountPath: '/workspace'] ]
    }

An equivalent declaration using the curly brackets notation is shown below::

    k8s {
        volumeClaims {
            'your-pvc-name' { mountPath = '/workspace' }
        }
    }

More than one volume claims can be defined repeating the name and mount path declaration in the ``volumeClaims`` block.


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
================== ================

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

.. _config-report:

Scope `report`
--------------

The ``report`` scope scope allows you to define configuration setting of the workflow :ref:`execution-report`.

================== ================
Name                Description
================== ================
enabled             If ``true`` it create the workflow execution report.
file                The path of the created execution report file (default: ``report.html``).
================== ================


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

.. tip:: Two or more configuration profiles can be specified by separating the profile names
    with a comma character, for example::

        nextflow run <your script> -profile standard,cloud

The above feature requires version 0.28.x or higher. 

Environment variables
=====================

The following environment variables control the configuration of the Nextflow runtime and
the Java virtual machine used by it.

=========================== ================
Name                        Description
=========================== ================
NXF_HOME                    Nextflow home directory (default: ``$HOME/.nextflow``).
NXF_VER                     Defines what version of Nextflow to use.
NXF_ORG                     Default `organization` prefix when looking for a hosted repository (default: ``nextflow-io``).
NXF_GRAB                    Provides extra runtime dependencies downloaded from a Maven repository service.
NXF_OPTS                    Provides extra options for the Java and Nextflow runtime. It must be a blank separated list of ``-Dkey[=value]`` properties.
NXF_CLASSPATH               Allows to extend the Java runtime classpath with extra jar files or class folders.
NXF_ASSETS                  Defined the directory where downloaded pipeline repositories are stored (default: ``$NXF_HOME/assets``)
NXF_PID_FILE                Name of the file where the process PID is saved when Nextflow is launched in background.
NXF_WORK                    Directory where working files are stored (usually your *scratch* directory)
NXF_TEMP                    Directory where temporary files are stored
NXF_DEBUG                   Defines scripts debugging level: ``1`` dump task environment variables in the task log file; ``2`` enables command script execution tracing; ``3`` enables command wrapper execution tracing.
NXF_EXECUTOR                Defines the default process executor e.g. `sge`
NXF_SINGULARITY_CACHEDIR    Directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible to all computing nodes.
NXF_JAVA_HOME               Defines the path location of the Java VM installation used to run Nextflow. This variable overrides the ``JAVA_HOME`` variable if defined.
JAVA_HOME                   Defines the path location of the Java VM installation used to run Nextflow.
JAVA_CMD                    Defines the path location of the Java binary command used to launch Nextflow.
HTTP_PROXY                  Defines the HTTP proxy server
HTTPS_PROXY                 Defines the HTTPS proxy server
=========================== ================
