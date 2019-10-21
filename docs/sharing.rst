.. _sharing-page:

****************
Pipeline sharing
****************

Nextflow seamlessly integrates with `BitBucket <http://bitbucket.org/>`_ [#]_, `GitHub <http://github.com>`_,
and `GitLab <http://gitlab.com>`_ hosted code repositories and sharing platforms. This feature allows you to manage your
project code in a more consistent manner or use other people's Nextflow pipelines, published through BitBucket/GitHub/GitLab,
in a quick and transparent way.

How it works
============

When you launch a script execution with Nextflow, it will look for a file with the pipeline name you've specified.
If that file does not exist, it will look for a public repository with the same name on GitHub (unless otherwise specified).
If it is found, the repository is automatically downloaded to your computer and executed. This repository is
stored in the Nextflow home directory, that is by default the ``$HOME/.nextflow`` path, and thus will be reused for any further
executions.

Running a pipeline
==================

To launch the execution of a pipeline project, hosted in a remote code repository, you simply need to specify its `qualified` name
or the repository URL after the ``run`` command. The qualified name is formed by two parts: the `owner` name and the
`repository` name separated by a ``/`` character.

In other words if a Nextflow project is hosted, for example, in a GitHub repository at the address
``http://github.com/foo/bar``, it can be executed by entering the following command in your shell terminal::

    nextflow run foo/bar

or using the project URL::

    nextflow run http://github.com/foo/bar


.. note:: In the first case, if your project is hosted on a service other than GitHub, you will need to specify this hosting
    service in the command line by using the ``-hub`` option. For example ``-hub bitbucket`` or ``-hub gitlab``.
    In the second case, i.e. when using the project URL as name, the ``-hub`` option is not needed.

You can try this feature out by simply entering the following command in your shell terminal::

    nextflow run nextflow-io/hello

It will download a trivial `Hello` example from the repository published at the following address
http://github.com/nextflow-io/hello and execute it in your computer.

If the `owner` part in the pipeline name is omitted, Nextflow will look for a pipeline between the ones you have
already executed having a name that matches the name specified. If none is found it will try to download
it using the `organisation` name defined by the environment variable ``NXF_ORG`` (which by default is ``nextflow-io``).


.. tip:: To access a private repository, specify the access credentials by using the ``-user`` command
    line option, then the program will ask you to enter the password interactively.
    Private repositories access credentials can also be defined in the `SCM configuration file`_.


Handling revisions
==================

Any Git branch, tag or commit ID defined in a project repository, can be used to specify the revision that you want to execute
when launching a pipeline by adding the ``-r`` option to the run command line. So for example you could enter::

    nextflow run nextflow-io/hello -r mybranch

or ::

    nextflow run nextflow-io/hello -r v1.1


It will execute two different project revisions corresponding to the Git tag/branch having that names.

Commands to manage projects
===========================

The following commands allows you to perform some basic operations that can be used to manage your projects.

.. note:: Nextflow is not meant to replace functionalities provided by the `Git <https://git-scm.com/>`_ tool. You may still need it to create new
  repositories or commit changes, etc.

Listing available projects
--------------------------

The ``list`` command allows you to list all the projects you have downloaded in your computer. For example::

    nextflow list

This prints a list similar to the following one::

    cbcrg/ampa-nf
    cbcrg/piper-nf
    nextflow-io/hello
    nextflow-io/examples


Showing project information
---------------------------

By using the ``info`` command you can show information from a downloaded project. For example::

     project name: nextflow-io/hello
     repository  : http://github.com/nextflow-io/hello
     local path  : $HOME/.nextflow/assets/nextflow-io/hello
     main script : main.nf
     revisions   :
     * master (default)
       mybranch
       v1.1 [t]
       v1.2 [t]

Starting from the top it shows: 1) the project name; 2) the Git repository URL; 3) the local folder where the
project has been downloaded; 4) the script that is executed when launched; 5) the list of available
revisions i.e. branches and tags. Tags are marked with a ``[t]`` on the right, the current checked-out revision is
marked with a ``*`` on the left.

Pulling or updating a project
-----------------------------

The ``pull`` command allows you to download a project from a GitHub repository or to update it if
that repository has already been downloaded. For example::

    nextflow pull nextflow-io/examples

Altenatively, you can use the repository URL as the name of the project to pull::

    nextflow pull https://github.com/nextflow-io/examples


Downloaded pipeline projects are stored in the folder ``$HOME/.nextflow/assets`` in your computer.


Viewing the project code
-------------------------

The ``view`` command allows you to quickly show the content of the pipeline script you have downloaded. For example::

    nextflow view nextflow-io/hello

By adding the ``-l`` option to the example above it will list the content of the repository.


Cloning a project into a folder
-------------------------------

The ``clone`` command allows you to copy a Nextflow pipeline project to a directory of your choice. For example::

    nextflow clone nextflow-io/hello target-dir

If the destination directory is omitted the specified project is cloned to a directory with the same name as the
pipeline base name (e.g. `hello`) in the current folder.

The clone command can be used to inspect or modify the source code of a pipeline project. You can eventually commit and push
back your changes by using the usual Git/GitHub workflow.

Deleting a downloaded project
-----------------------------

Downloaded pipelines can be deleted by using the ``drop`` command, as shown below::

    nextflow drop nextflow-io/hello

SCM configuration file
=======================

The file ``$HOME/.nextflow/scm`` allows you to centralise the security credentials required to access private project
repositories on Bitbucket, GitHub and GitLab source code management (`SCM`) platforms or to manage the configuration properties
of private server installations (of the same platforms).

The configuration properties for each SCM platform are defined inside the ``providers`` section,
properties for the same provider are grouped together with a common name and delimited with curly brackets as in this example::

    providers {
        <provider-name> {
            property = value
            :
        }
    }


In the above template replace `<provider-name>` with one of the "default" servers (i.e. ``bitbucket``, ``github`` or ``gitlab``)
or a custom identifier representing a private SCM server installation.

The following configuration properties are supported for each provider configuration:

=================== ==============
Name                Description
=================== ==============
user                User name required to access private repositories on the SCM server.
password            User password required to access private repositories on the SCM server.
token               Private API access token (used only when the specified platform is ``gitlab``).
:sup:`*` platform   SCM platform name, either: ``github``, ``gitlab`` or ``bitbucket``.
:sup:`*` server     SCM server name including the protocol prefix e.g. ``https://github.com``.
:sup:`*` endpoint   SCM API `endpoint` URL e.g. ``https://api.github.com`` (default: the same value specified for ``server``).
=================== ==============

The attributes marked with a * are only required when defining the configuration of a private SCM server.


BitBucket credentials
---------------------

Create a ``bitbucket`` entry in the `SCM configuration file`_ specifying your user name and password, as shown below::

    providers {

        bitbucket {
            user = 'me'
            password = 'my-secret'
        }

    }

BitBucket Server credentials
-----------------------------

`BitBucket Server <https://confluence.atlassian.com/bitbucketserver>`_ is a self-hosted Git repository and management
platform.

.. note::
    BitBucket Server uses different API from the `BitBucket <https://bitbucket.org/>`_ cloud service. Make sure to
    use the right configuration whether you are using the cloud service or a self-hosted installation.

To access your local BitBucket Server create an entry in the `SCM configuration file`_ specifying as shown below::

        providers {

            mybitbucket {
                platform = 'bitbucketserver'
                server = 'https://your.bitbucket.host.com'
                endpoint = 'https://your.bitbucket.host.com'
                user = 'your-user'
                password = 'your-password or your-token'
            }

        }


GitHub credentials
------------------

Create a ``github`` entry in the `SCM configuration file`_ specifying your user name and password as shown below::

    providers {

        github {
            user = 'me'
            password = 'my-secret'
        }

    }

.. tip:: You can use use a `Personal API token <https://github.com/blog/1509-personal-api-tokens>`_ in place of your
    GitHub password.


GitLab credentials
-------------------

Create a ``gitlab`` entry in the `SCM configuration file`_ specifying the user name, password and your API access token
that can be found in your GitLab `account page <https://gitlab.com/profile/account>`_ (sign in required). For example::

    providers {

        gitlab {
            user = 'me'
            password = 'my-secret'
            token = 'YgpR8m7viH_ZYnC8YSe8'
        }

    }


Gitea credentials
-----------------

`Gitea <https://gitea.io>`_ is a Git repository server with GitHub-like GUI access. Since Gitea installation is quite 
easy, it is suitable for building a private development environment in your network. To access your Gitea server, you 
have to provide all the credential information below::

    providers {

        mygitea {
            server = 'http://your-domain.org/gitea'
            endpoint = 'http://your-domain.org/gitea/api/v1'
            platform = 'gitea'
            user = 'your-user'
            password = 'your-password'
            token = 'your-api-token'
        }

    }


See `Gitea documentation <https://docs.gitea.io/en-us/api-usage/>`_ about how to enable API access on your 
server and how to issue a token. 


Private server configuration
============================

Nextflow is able to access repositories hosted on private BitBucket, GitHub, GitLab and Gitea server installations.

In order to use a private SCM installation you will need to set the server name and access credentials
in your `SCM configuration file`_ .

If, for example, the host name of your private GitLab server is ``gitlab.acme.org``, you will need to have in the
``$HOME/.nextflow/scm`` file a configuration like the following::

    providers {

        mygit {
            server = 'http://gitlab.acme.org'
            platform = 'gitlab'
            user = 'your-user'
            password = 'your-password'
            token = 'your-api-token'
        }

    }


Then you will be able to run/pull a project with Nextflow using the following command line::

    $ nextflow run foo/bar -hub mygit

Or, in alternative, using the Git clone URL::

    $ nextflow run http://gitlab.acme.org/foo/bar.git

.. warning:: When accessing a private SCM installation over ``https`` and that server uses a custom SSL certificate
  you may need to import such certificate into your local Java keystore. Read more
  `here <https://docs.oracle.com/javase/tutorial/security/toolsign/rstep2.html>`_.

Local repository configuration
==============================

Nextflow is also able to handle repositories stored in a local or shared file system. The repository
must be created as a `bare repository <https://mijingo.com/blog/what-is-a-bare-git-repository>`_.


Having, for example. a bare repository store at path ``/shared/projects/foo.git``, Nextflow is able
to run it using the following syntax::

  $ nextflow run file:/shared/projects/foo.git

See `Git documentation <https://git-scm.com/book/en/v2/Git-on-the-Server-Getting-Git-on-a-Server>`_ for
more details about how create and manage bare repositories.

Publishing your pipeline
========================

In order to publish your Nextflow pipeline to GitHub (or any other supported platform) and allow other people to use it,
you only need to create a GitHub repository containing all your project script and data files. If you don't know how to
do it, follow this simple tutorial that explains how `create a GitHub repository <https://help.github.com/articles/create-a-repo>`_.

Nextflow only requires that the main script in your pipeline project is called ``main.nf``. A different name can be
used by specifying the ``manifest.mainScript`` attribute in the ``nextflow.config`` file that must be
included in your project. For example::

  manifest.mainScript = 'my_very_long_script_name.nf'

To learn more about this and other project metadata information, that can be defined in the Nextflow configuration file,
read the :ref:`Manifest <config-manifest>` section on the Nextflow configuration page.

Once you have uploaded your pipeline project to GitHub other people can execute it simply
using the project name or the repository URL.

For if your GitHub account name is ``foo`` and you have uploaded a project into a repository named ``bar`` the
repository URL will be ``http://github.com/foo/bar`` and people will able to download and run it by using either
the command::

    nextflow run foo/bar

or

::

    nextflow run http://github.com/foo/bar

See the `Running a pipeline`_ section for more details on how to run Nextflow projects.

Manage dependencies
=====================

Computational pipelines are rarely composed by a single script. In real world applications they depend on dozens of other components.
These can be other scripts, databases, or applications compiled for a platform native binary format.

External dependencies are the most common source of problems when sharing a piece of software, because the
users need to have an identical set of tools and the same configuration to be able to use it. In many cases this has proven to be
a painful and error prone process, that can severely limit the ability to reproduce computational results on a system other than
the one on which it was originally developed.

Nextflow tackles this problem by integrating GitHub, BitBucket and GitLab sharing platforms and
`Docker <http://www.docker.com>`_ containers technology.

The use of a code management system is important to keep together all the dependencies of your
pipeline project and allows you to track the changes of the source code in a consistent manner.

Moreover to guarantee that a pipeline is reproducible it should be self-contained i.e. it should have ideally no
dependencies on the hosting environment. By using Nextflow you can achieve this goal following these methods:

Third party scripts
--------------------

Any third party script that does not need to be compiled (Bash, Python, Perl, etc) can be included in the pipeline
project repository, so that they are distributed with it.

Grant the execute permission to these files and copy them into a folder named ``bin/`` in the root directory of your
project repository. Nextflow will automatically add this folder to the ``PATH`` environment variable, and the scripts
will automatically be accessible in your pipeline without the need to specify an absolute path to invoke them.

System environment
--------------------

Any environment variable that may be required by the tools in your pipeline can be defined in the ``nextflow.config`` file
by using the ``env`` scope and including it in the root directory of your project. For example::

  env {
    DELTA = 'foo'
    GAMMA = 'bar'
  }


See the :ref:`config-page` page to learn more about the Nextflow configuration file.

Resource manager
--------------------

When using Nextflow you don't need to write the code to parallelize your pipeline for a specific grid engine/resource
manager because the parallelization is defined implicitly and managed by the Nextflow runtime. The target execution
environment is parametrized and defined in the configuration file, thus your code is free from this kind of dependency.

Bootstrap data
--------------------

Whenever your pipeline requires some files or dataset to carry out any initialization step, you
can include this data in the pipeline repository itself and distribute them together.

To reference this data in your pipeline script in a portable manner (i.e. without the need to use a static absolute path)
use the implicit variable ``baseDir`` which locates the base directory of your pipeline project.

For example, you can create a folder named ``dataset/`` in your repository root directory and copy there the
required data file(s) you may need, then you can access this data in your script by writing::

   sequences = file("$baseDir/dataset/sequences.fa")
   sequences.splitFasta {
        println it
    }

User inputs
-------------

Nextflow scripts can be easily parametrised to allow users to provide their own input data. Simply declare on the
top of your script all the parameters it may require as shown below::

  params.my_input = 'default input file'
  params.my_output = 'default output path'
  params.my_flag = false
  ..

The actual parameter values can be provided when launching the script execution on the command line
by prefixed the parameter name with a double minus character i.e. ``--``, for example::

  nextflow run <your pipeline> --my_input /path/to/input/file --my_output /other/path --my_flag true




Binary applications
--------------------

Docker allows you to ship any binary dependencies that you may have in your pipeline to a portable image
that is downloaded on-demand and can be executed on any platform where a Docker engine is installed.

In order to use it with Nextflow, create a Docker image containing the tools needed by your pipeline and make it available
in the `Docker registry <https://registry.hub.docker.com>`_.

Then declare in the ``nextflow.config`` file, that you will include in your project, the name of the Docker image you
have created. For example::

  process.container = 'my-docker-image'
  docker.enabled = true

In this way when you launch the pipeline execution, the Docker image will be automatically downloaded and used to run 
your tasks.

Read the :ref:`docker-page` page to lean more on how to use Docker containers with Nextflow.


This mix of technologies makes it possible to write self-contained and truly reproducible pipelines which require
zero configuration and can be reproduced in any system having a Java VM and a Docker engine installed.


.. [#] BitBucket provides two types of version control system: `Git` and `Mercurial`. Nextflow supports only `Git` based repositories.
