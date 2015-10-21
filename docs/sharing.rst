.. _sharing-page:

****************************
Pipeline sharing
****************************

Nextflow seamlessly integrates with `GitHub <http://github.com>`_, `BitBucket <http://bitbucket.org/>`_ [#]_
and `GitLab <http://gitlab.com>`_ hosted code repositories and sharing platforms. This feature allows you to manage your
project code in a more consistent manner, or use other people's Nextflow pipelines, published through GitHub/BitBucket/GitLab,
in a quick and transparent way.

How it works
=============

When you launch a script execution with Nextflow, it will look for a file with the pipeline name you've specified.
If that file does not exist, it will look for a public repository with the same name on GitHub (unless otherwise specified).
If it is found, the repository is automatically downloaded to your computer and the code executed. This repository is
stored in the Nextflow home directory, by default the ``$HOME/.nextflow`` path, thus it will be reused for any further
execution.

Run a pipeline
================

To launch the execution of a pipeline project hosted in remote code repository you simply need to specify its `qualified` name
after the ``run`` command or the repository URL. The qualified name is formed by two parts: the `owner` name and the `repository` name separated by a ``/`` character.

In other words if a Nextflow project is hosted, for example, into a GitHub repository having an address like
``http://github.com/foo/bar``, It can be executed by entering in your shell terminal the following command::

    nextflow run foo/bar

or using the project URL::

    nextflow run http://github.com/foo/bar


.. note:: In the first case, if your project is hosted on a service other than GitHub, you will need to specify this hosting
    service in the command line by using the ``-hub`` command line option. For example ``-hub bitbucket`` or ``-hub gitlab``.
    In the second case, i.e., when using the project URL as name, that is not needed.

.. tip:: Private repositories can be accessed by specifying the access credentials by using the ``-user`` command
  line option, then the program will ask to enter the password interactively.


You can try this feature out, by simply entering the following command in your shell terminal::

    nextflow run nextflow-io/hello

It will download a trivial `Hello` example from the repository published at the following address
http://github.com/nextflow-io/hello and execute it in your computer.

If `owner` part in the pipeline name is omitted, Nextflow will look for a pipeline between the ones you have
already executed having a name that matches the name specified. If none it's found it will try to download
it using the `organisation` name defined by the environment variable ``NXF_ORG`` (which by default is ``nextflow-io``).

Handle revisions
==================

Any branch, tag or commit ID defined in a project Git repository can be used to specify a revision, that you want to execute,
when running your pipeline by adding the ``-r`` option to the run command line. So for example you could enter::

    nextflow run nextflow-io/hello -r mybranch

or ::

    nextflow run nextflow-io/hello -r v1.1


It will execute two different code revisions as specified.

Commands to manage projects
============================

The following commands allows you to perform some basic operations that can be used to manage your projects.

.. note:: Anyway Nextflow is not meant to replace functionalities provided by the `Git` tool, you may still need it to create new
  repositories or commit changes, etc.

List available projects
-------------------------

The ``ls`` command allows you to list all the pipeline projects you have downloaded in your computer. For example::

    nextflow ls

This prints a list similar to the following one::

    cbcrg/piper-nf
    nextflow-io/hello


Show project information
--------------------------

By using the ``info`` command you can show information from a downloaded pipeline. For example::

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
pipeline has been downloaded; 4) the script that is executed when launched; 5) the list of available
revisions i.e. branches and tags. Tags are marked with a ``[t]`` on the right, the current checked-out revision is
marked with a ``*`` on the left.

Pull or update a project
--------------------------

The ``pull`` command allows you to download a project from a GitHub repository or to update it if
that repository has already been downloaded. For example::

    nextflow pull nextflow-io/examples

Downloaded pipeline projects are stored in the folder ``$HOME/.nextflow/assets`` in your computer.


View the project code
--------------------------

The ``view`` command allows you to show quickly the content of a pipeline script you may have downloaded. For example::

    nextflow view nextflow-io/hello

By adding the ``-l`` option to the example above it will list the content of the repository.


Clone a project into a folder
-------------------------------

The ``clone`` command allows you to copy a Nextflow pipeline project to a directory of your choice. For example::

    nextflow clone nextflow-io/hello target-dir

If the destination directory is omitted the specified project is cloned to a directory with the same name as the
pipeline base name (e.g. hello) in the current folder.

The clone command can be used to inspect or modify the source code of a pipeline project. You can eventually commit and push
back your changes by using the usual Git/GitHub workflow.

Delete a downloaded project
-----------------------------

Downloaded pipelines can be deleted by using the ``drop`` command, as shown below::

    nextflow drop nextflow-io/hello

SCM configuration file
=======================

The file ``$HOME/.nextflow/scm`` allows you to centralise the security credentials required to access private project
repositories on GitHub, Bitbucket and GitLab source code management (`SCM`) platforms or to manage the configuration properties
of private server installations (of the same platforms).

The configuration properties for each SCM platform are defined inside the ``providers`` section,
properties for the same provider are grouped together with a common name and delimited with curly brackets as in this example::

    providers {
        <provider-name> {
            property = value
            :
        }
    }


Where `<provider-name>` need to be replaced with one of the "default" server (i.e. ``github``, ``gitlab`` and ``bitbucket``)
or a custom identifier representing a private SCM server installation.

The following configuration properties are supported for each provider configuration:

=================== ==============
Name                Description
=================== ==============
user                User name required to access private repositories on the SCM server.
password            User password required to access private repository on the SCM server.
token               Private API access token (used only when the specified platform is ``gitlab``).
:sup:`*` platform   SCM platform name, either: ``github``, ``gitlab`` or ``bitbucket``.
:sup:`*` server     SCM server name including the protocol prefix e.g. ``https://github.com``.
:sup:`*` endpoint   SCM API `endpoint` URL e.g. ``https://api.github.com`` (default: the same value provided for ``server``).
=================== ==============

The attribute marked with a * are only required when defining the configuration of a private SCM servers.


GitHub credentials
-------------------

Create an ``github`` entry in the `SCM configuration file`_ specifying your user name and password as shown below::

    providers {

        github {
            user = 'me'
            password = 'my-secret'
        }

    }

.. tip:: You can use use a `Personal API token <https://github.com/blog/1509-personal-api-tokens>`_ in place of your GitHub password.


Bitbucket credentials
----------------------

Create an ``bitbucket`` entry in the `SCM configuration file`_ specifying your the user name and password, as shown below::

    providers {

        bitbucket {
            user = 'me'
            password = 'my-secret'
        }

    }


GitLab credentials
-------------------

Create an ``bitbucket`` entry in the `SCM configuration file`_ specifying the user name, password and your API access token
that can be found in your GitLab `account page <https://gitlab.com/profile/account>`_ (sign in required). For example::

    providers {

        gitlab {
            user = 'me'
            password = 'my-secret'
            token = 'YgpR8m7viH_ZYnC8YSe8'
        }

    }


Private server configuration
=============================

Nextflow is able to access repositories hosted on private server installations of the GitHub, Bitbucket and GitLab SCM platforms.

In order of being able to use a private SCM server you will need to include the server name and access credentials of your SCM installation
in your `SCM configuration file`_ .

If, for example, the host name of your private GitLab server is ``github.acme.org``, you will need to to add in the
``$HOME/.nextflow/scm`` file a configuration like the following one::

    providers {

        mygit {
            server = 'http://github.acme.org'
            platform = 'github'
            user = 'your-user'
            password = 'your-password'
        }

    }


Then you will be able to run/pull project with Nextflow using the following command line::

    $ nextflow run foo/bar -hub mygit

Or, in alternative, using the Git clone URL::

    $ nextflow run http://github.acme.org/foo/bar.git


Publish your pipeline
======================

In order to publish your Nextflow pipeline to GitHub and allow other people to use it, you simply need to create a
GitHub repository containing all your project script and data files. If you don't know how to do it, follow this simple tutorial
that explains how `create a GitHub repository <https://help.github.com/articles/create-a-repo>`_.

Nextflow only requires that main script in your pipeline project is called ``main.nf``. A different name can be
used by specifying the ``manifest.mainScript`` attribute in the ``nextflow.config`` file that must be
included to your project. For example::

  manifest.mainScript = 'my_pipeline_very_long_name.nf'

Learn more about this and other project meta-data information that can be defined in the Nextflow configuration file
read the :ref:`Manifest <config-manifest>` section in the Nextflow configuration page.

Once you have uploaded your pipeline project into GitHub (or BitBucket/GitLab) other people can use it by specifying the
pipeline `qualified` name on the Nextflow run command line. The qualified name is simply the GitHub user name
(or organisation) plus the repository name.

For the sake of the example if your GitHub account name is ``foo`` and you have uploaded it to a repository named ``bar`` the
repository URL will be ``http://github.com/foo/bar`` and people will able to run it by entering the command::

  nextflow run foo/bar



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

Any third party script that does not need to be compiled (BASH, Python, Perl, etc) can be included in the pipeline
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
