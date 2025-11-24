(sharing-page)=

# Sharing pipelines

<<<<<<< HEAD
Nextflow seamlessly integrates with [BitBucket](http://bitbucket.org/) [^id2], [GitHub](http://github.com), and [GitLab](http://gitlab.com) hosted code repositories and sharing platforms. This feature allows you to manage your project code in a more consistent manner or use other people's Nextflow pipelines, published through BitBucket/GitHub/GitLab, in a quick and transparent way.

## How it works

When you launch a script execution with Nextflow, it will look for a file with the pipeline name you've specified. If that file does not exist, it will look for a public repository with the same name on GitHub (unless otherwise specified). If it is found, the repository is automatically downloaded to your computer and executed. This repository is stored in the Nextflow home directory, that is by default the `$HOME/.nextflow` path, and thus will be reused for any further executions.

## Running a pipeline

To launch the execution of a pipeline project, hosted in a remote code repository, you simply need to specify its qualified name or the repository URL after the `run` command. The qualified name is formed by two parts: the `owner` name and the `repository` name separated by a `/` character.

In other words if a Nextflow project is hosted, for example, in a GitHub repository at the address `http://github.com/foo/bar`, it can be executed by entering the following command in your shell terminal:

```bash
nextflow run foo/bar
```

or using the project URL:

```bash
nextflow run http://github.com/foo/bar
```

:::{note}
In the first case, if your project is hosted on a service other than GitHub, you will need to specify this hosting service in the command line by using the `-hub` option. For example `-hub bitbucket` or `-hub gitlab`. In the second case, i.e. when using the project URL as name, the `-hub` option is not needed.
:::

You can try this feature out by simply entering the following command in your shell terminal:

```bash
nextflow run nextflow-io/hello
```

It will download a trivial `Hello` example from the repository published at the following address <http://github.com/nextflow-io/hello> and execute it in your computer.

If the `owner` part in the pipeline name is omitted, Nextflow will look for a pipeline between the ones you have already executed having a name that matches the name specified. If none is found it will try to download it using the `organisation` name defined by the environment variable `NXF_ORG` (which by default is `nextflow-io`).

:::{tip}
To access a private repository, specify the access credentials by using the `-user` command line option, then the program will ask you to enter the password interactively. Private repository access credentials can also be defined in the [SCM configuration file](#scm-configuration-file)(#s.
:::

## Handling revisions

Any Git branch, tag or commit ID defined in a project repository, can be used to specify the revision that you want to execute when launching a pipeline by adding the `-r` option to the run command line. So for example you could enter:

```bash
nextflow run nextflow-io/hello -r mybranch
```

or

```bash
nextflow run nextflow-io/hello -r v1.1
```

It will execute two different project revisions corresponding to the Git tag/branch having that names.

:::{versionadded} 24.XX.0-edge
:::

Nextflow downloads and locally maintains each explicitly requested Git branch, tag or commit ID in a separate directory path, thus enabling to run multiple revisions of the same pipeline at the same time. Each downloaded revision is stored in a subdirecrory path of the local project path: `$NXF_ASSETS/<org>/<repo>/.nextflow/commits/<commitId>`.

:::{warning}
If you really care about reproducibility of your pipelines, consider explicitly referring to them by tag or commit ID, rather than my branch. This is because the same branch will point to different underlying commits over time, as pipeline development goes on.
:::

## Commands to manage projects

The following commands allows you to perform some basic operations that can be used to manage your projects.
=======
Nextflow seamlessly integrates with popular Git providers, including [BitBucket](http://bitbucket.org/), [GitHub](http://github.com), and [GitLab](http://gitlab.com) for managing Nextflow pipelines as version-controlled Git repositories. This feature allows you to easily use other people's Nextflow pipelines and publish your own pipelines.
>>>>>>> master

:::{note}
Nextflow is not meant to completely replace the [Git](https://git-scm.com/) tool. You may still need `git` to create new repositories or commit changes, etc.
:::

## Git configuration

You can configure your credentials for various Git providers in the Git configuration file, located at `$HOME/.nextflow/scm`. See {ref}`git-page` for more information.

## Using a local repository

Nextflow can work with repositories stored in a local or shared file system. The repository must be created as a [bare repository](https://craftquest.io/articles/what-is-a-bare-git-repository).

For example, given a bare repository at `/shared/projects/hello.git`, Nextflow is able to run it using the following syntax:

```bash
<<<<<<< HEAD
nextflow list
```

This prints a list similar to the following one:

```
cbcrg/ampa-nf
cbcrg/piper-nf
nextflow-io/hello
nextflow-io/examples
```

### Showing project information

By using the `info` command you can show information from a downloaded project. For example:

```console
$ nextflow info hello
project name: nextflow-io/hello
repository  : http://github.com/nextflow-io/hello
local path  : /Users/evanfloden/.nextflow/assets/nextflow-io/hello
main script : main.nf
revisions   :
P master (default)
  mybranch
P v1.1 [t]
  v1.2 [t]
```

Starting from the top it shows: 1) the project name; 2) the Git repository URL; 3) the local path where the project can be found; 4) the script that is executed when launched; 5) the list of available revisions i.e. branches and tags. Tags are marked with a `[t]` on the right, the locally pulled revisions are marked with a `P` on the left.

### Pulling or updating a project

The `pull` command allows you to download a project from a GitHub repository or to update it if that repository has already been downloaded. For example:

```bash
nextflow pull nextflow-io/examples
```

Alternatively, you can use the repository URL as the name of the project to pull:

```bash
nextflow pull https://github.com/nextflow-io/examples
```

Downloaded pipeline projects are stored in the folder `$HOME/.nextflow/assets` in your computer.

### Viewing the project code

The `view` command allows you to quickly show the content of the pipeline script you have downloaded. For example:

```bash
nextflow view nextflow-io/hello
```

By adding the `-l` option to the example above it will list the content of the repository.

### Cloning a project into a folder

The `clone` command allows you to copy a Nextflow pipeline project to a directory of your choice. For example:

```bash
nextflow clone nextflow-io/hello target-dir
```

If the destination directory is omitted the specified project is cloned to a directory with the same name as the pipeline base name (e.g. `hello`) in the current folder.

The clone command can be used to inspect or modify the source code of a pipeline project. You can eventually commit and push back your changes by using the usual Git/GitHub workflow.

### Deleting a downloaded project

Downloaded pipelines can be deleted by using the `drop` command, as shown below:

```bash
nextflow drop nextflow-io/hello
```

(sharing-scm-file)=

## SCM configuration file

The file `$HOME/.nextflow/scm` allows you to centralise the security credentials required to access private project repositories on Bitbucket, GitHub and GitLab source code management (SCM) platforms or to manage the configuration properties of private server installations (of the same platforms).

The configuration properties for each SCM platform are defined inside the `providers` section, properties for the same provider are grouped together with a common name and delimited with curly brackets as in this example:

```groovy
providers {
    <provider-name> {
        property = value
        // ...
    }
}
```

In the above template replace `<provider-name>` with one of the "default" servers (i.e. `bitbucket`, `github` or `gitlab`) or a custom identifier representing a private SCM server installation.

:::{versionadded} 20.10.0
A custom location for the SCM file can be specified using the `NXF_SCM_FILE` environment variable.
:::

The following configuration properties are supported for each provider configuration:

`providers.<provider>.user`
: User name required to access private repositories on the SCM server.

`providers.<provider>.password`
: User password required to access private repositories on the SCM server.

`providers.<provider>.token`
: *Required only for private Gitlab servers*
: Private API access token.

`providers.<provider>.platform`
: *Required only for private SCM servers*
: SCM platform name, either: `github`, `gitlab` or `bitbucket`.

`providers.<provider>.server`
: *Required only for private SCM servers*
: SCM server name including the protocol prefix e.g. `https://github.com`.

`providers.<provider>.endpoint`
: *Required only for private SCM servers*
: SCM API `endpoint` URL e.g. `https://api.github.com` (default: the same as `providers.<provider>.server`).

## SCM providers

### BitBucket credentials

Create a `bitbucket` entry in the [SCM configuration file](#scm-configuration-file) specifying your user name and app password, as shown below:

```groovy
providers {
    bitbucket {
        user = 'me'
        password = 'my-secret'
    }
}
```

:::{note}
App passwords are substitute passwords for a user account which you can use for scripts and integrating tools in order to avoid putting your real password into configuration files. Learn more at [this link](https://support.atlassian.com/bitbucket-cloud/docs/app-passwords/).
:::

### BitBucket Server credentials

[BitBucket Server](https://confluence.atlassian.com/bitbucketserver) is a self-hosted Git repository and management platform.

:::{note}
BitBucket Server uses a different API from the [BitBucket](https://bitbucket.org/) cloud service. Make sure to use the right configuration whether you are using the cloud service or a self-hosted installation.
:::

To access your local BitBucket Server create an entry in the [SCM configuration file](#scm-configuration-file) specifying as shown below:

```groovy
providers {
    mybitbucket {
        platform = 'bitbucketserver'
        server = 'https://your.bitbucket.host.com'
        endpoint = 'https://your.bitbucket.host.com'
        user = 'your-user'
        password = 'your-password or your-token'
    }
}
```

### GitHub credentials

Create a `github` entry in the [SCM configuration file](#scm-configuration-file) specifying your user name and access token as shown below:

```groovy
providers {
    github {
        user = 'your-user-name'
        password = 'your-personal-access-token'
    }
}
```

GitHub requires the use of a personal access token (PAT) in place of a password when accessing APIs. Learn more about PAT and how to create it at [this link](https://docs.github.com/en/github/authenticating-to-github/keeping-your-account-and-data-secure/creating-a-personal-access-token).

:::{versionadded} 23.01.0-edge
Nextflow automatically uses the `GITHUB_TOKEN` environment variable to authenticate access to the GitHub repository if no credentials are provided via the `scm` file. This is useful especially when accessing pipeline code from a GitHub Action. Read more about the token authentication in the [GitHub documentation](https://docs.github.com/en/actions/security-guides/automatic-token-authentication).
:::

### GitLab credentials

Create a `gitlab` entry in the [SCM configuration file](#scm-configuration-file) specifying the user name, password and your API access token that can be found in your GitLab [account page](https://gitlab.com/profile/account) (sign in required). For example:

```groovy
providers {
    gitlab {
        user = 'me'
        password = 'my-secret'
        token = 'YgpR8m7viH_ZYnC8YSe8'
    }
}
```

:::{tip}
The GitLab *token* string can be used as the `password` value in the above setting. When doing that the `token` field can be omitted.
:::

### Gitea credentials

[Gitea](https://gitea.io) is a Git repository server with GitHub-like GUI access. Since Gitea installation is quite easy, it is suitable for building a private development environment in your network. To access your Gitea server, you have to provide all the credential information below:

```groovy
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
```

See [Gitea documentation](https://docs.gitea.io/en-us/api-usage/) about how to enable API access on your server and how to issue a token.

### Azure Repos credentials

Nextflow has a builtin support for [Azure Repos](https://azure.microsoft.com/en-us/services/devops/repos/), a Git source code management service hosted in the Azure cloud. To access your Azure Repos with Nextflow provide the repository credentials using the configuration snippet shown below:

```groovy
providers {
    azurerepos {
        user = 'your-user-name'
        password = 'your-personal-access-token'
    }
}
```

:::{tip}
The Personal access token can be generated in the repository `Clone Repository` dialog.
:::

(aws-codecommit)=

### AWS CodeCommit credentials

:::{versionadded} 22.06.0-edge
:::

Nextflow supports [AWS CodeCommit](https://aws.amazon.com/codecommit/) as a Git provider to access and to share pipelines code.

To access your project hosted on AWS CodeCommit with Nextflow provide the repository credentials using the configuration snippet shown below:

```groovy
providers {
    my_aws_repo {
        platform = 'codecommit'
        user = '<AWS ACCESS KEY>'
        password = '<AWS SECRET KEY>'
    }
}
```

In the above snippet replace `<AWS ACCESS KEY>` and `<AWS SECRET KEY>` with your AWS credentials, and `my_aws_repo` with a name of your choice.

:::{tip}
The `user` and `password` are optional settings, if omitted the [AWS default credentials provider chain](https://docs.aws.amazon.com/sdk-for-java/v1/developer-guide/credentials.html) is used.
:::

Then the pipeline can be accessed with Nextflow as shown below:

```bash
nextflow run https://git-codecommit.eu-west-1.amazonaws.com/v1/repos/my-repo
```

In the above example replace `my-repo` with your own repository. Note also that AWS CodeCommit has different URLs depending the region in which you are working.

:::{note}
The support for protocols other than HTTPS is not available at this time.
:::

## Private server configuration

Nextflow is able to access repositories hosted on private BitBucket, GitHub, GitLab and Gitea server installations.

In order to use a private SCM installation you will need to set the server name and access credentials in your [SCM configuration file](#scm-configuration-file) .

If, for example, the host name of your private GitLab server is `gitlab.acme.org`, you will need to have in the `$HOME/.nextflow/scm` file a configuration like the following:

```groovy
providers {
    mygit {
        server = 'http://gitlab.acme.org'
        platform = 'gitlab'
        user = 'your-user'
        password = 'your-password'
        token = 'your-api-token'
    }
}
```

Then you will be able to run/pull a project with Nextflow using the following command line:

```bash
nextflow run foo/bar -hub mygit
```

Or, in alternative, using the Git clone URL:

```bash
nextflow run http://gitlab.acme.org/foo/bar.git
```

:::{note}
You must also specify the server API endpoint URL if it differs from the server base URL. For example, for GitHub Enterprise V3, add `endpoint = 'https://git.your-domain.com/api/v3'`.
:::

:::{warning}
When accessing a private SCM installation over `https` from a server that uses a custom SSL certificate, you may need to import the certificate into your local Java keystore. Read more [here](https://docs.oracle.com/javase/tutorial/security/toolsign/rstep2.html).
:::

## Local repository configuration

Nextflow is also able to handle repositories stored in a local or shared file system. The repository must be created as a [bare repository](https://mijingo.com/blog/what-is-a-bare-git-repository).

Having, for example. a bare repository store at path `/shared/projects/foo.git`, Nextflow is able to run it using the following syntax:

```bash
nextflow run file:/shared/projects/foo.git
=======
nextflow run file:/shared/projects/hello.git
>>>>>>> master
```

See [Git documentation](https://git-scm.com/book/en/v2/Git-on-the-Server-Getting-Git-on-a-Server) for more details about how create and manage bare repositories.

## Publishing your pipeline

In order to publish your Nextflow pipeline to GitHub (or any other supported platform) and allow other people to use it, you only need to create a GitHub repository containing all your project script and data files. If you don't know how to do it, follow this simple tutorial that explains how [create a GitHub repository](https://help.github.com/articles/create-a-repo).

Nextflow only requires that the main script in your pipeline project is called `main.nf`. A different name can be used by specifying the `manifest.mainScript` attribute in the `nextflow.config` file that must be included in your project. For example:

```groovy
manifest.mainScript = 'my_very_long_script_name.nf'
```

To learn more about this and other project metadata information, that can be defined in the Nextflow configuration file, read the {ref}`Manifest <config-manifest>` section on the Nextflow configuration page.

Once you have uploaded your pipeline project to GitHub other people can execute it simply using the project name or the repository URL.

For if your GitHub account name is `acme` and you have uploaded a project into a repository named `hello` the repository URL will be `http://github.com/acme/hello` and people will able to download and run it by using either the command:

```bash
nextflow run acme/hello
```

or

```bash
nextflow run http://github.com/acme/hello
```

See the {ref}`CLI <cli-page>` page to learn how to use the Nextflow command line to run pipelines and manage pipeline projects.

## Managing dependencies

Computational pipelines are rarely composed by a single script. In real world applications they depend on many other components, including other scripts and tools, databases, and specialized environments which provide compute and storage.

These external dependencies are the primary challenge when sharing software, because the users need to recreate the environment around a tool in order to use it. This setup process is often painful and error prone, which severely hinders the ability to reproduce computational results on a system other than the one on which it was originally developed.

Nextflow tackles this problem by integrating with existing tools for reproducible software, namely Git for source code and Docker for containers. These tools allow you to keep all the dependencies of your pipeline project in one place and track changes over time with version control.

By making your pipeline project is self-contained, meaning all of its dependencies are fully defined in the project itself, you gain two major advantages:

- **Portability**: the pipeline can be run in virtually any environment with a Java VM and a container runtime
- **Reproducibility**: any results produced by the pipelined can be easily reproduced, even across different environments

One way to account for dependencies is to break them down into three categories: code, data, and environment. Here we will describe how to include each of these dependencies in your Nextflow pipeline:

### Code

Aside from pipeline scripts, you may have additional scripts and tools used by individual tasks.

#### Standard dependencies

Many standard tools can be accessed as Docker containers or Conda/Spack packages. In this case, you only need to specify the container image URL (e.g. from [DockerHub](https://hub.docker.com)) or package name with the `container` or `conda` directive, respectively.

Make sure to enable the desired method in your Nextflow configuration:

```groovy
// containers
docker.enabled = true

// conda packages
conda.enabled = true
```

This way, when you launch your pipeline, Nextflow will automatically download the necessary dependencies to run your tasks based on this configuration.

Read the {ref}`container-page` page to learn more about how to use containers with Nextflow, and the {ref}`conda-page` page for Conda packages.

:::{tip}
For maximal reproducibility, make sure to define a specific version for each tool. Otherwise, your pipeline might use different versions across subsequent runs, which can introduce subtle differences to your results.
:::

(bundling-executables)=

#### The `bin` directory

As for custom scripts, you can include executable scripts in the `bin` directory of your pipeline repository. When configured correctly, these scripts can be executed like a regular command from any process script (i.e. without modifying the `PATH` environment variable or using an absolute path), and changing the script will cause the task to be re-executed on a resumed run (i.e. just like changing the process script itself).

To configure a custom script:

1. Save the script in the `bin` directory (relative to the pipeline repository root).
2. Specify a portable shebang (see note below for details).
3. Make the script executable. For example: `chmod a+x bin/my_script.py`

:::{tip}
To maximize the portability of your bundled script, use `env` to dynamically resolve the location of the interpreter instead of hard-coding it in the shebang line.

For example, shebang definitions `#!/usr/bin/python` and `#!/usr/local/bin/python` both hard-code specific paths to the Python interpreter. Instead, the following approach is more portable:

```bash
#!/usr/bin/env python
```
:::

(lib-directory)=

#### The `lib` directory

The `lib` directory can be used to add utility code or external libraries without cluttering the pipeline scripts. The `lib` directory in the Nextflow project root is added to the classpath by default. Classes defined in the `lib` directory will be available in pipeline scripts. Functions defined outside of classes will not be available in pipeline scripts.

For example, `lib/DNASequence.groovy` defines the `DNASequence` class:

```groovy
// lib/DNASequence.groovy
class DNASequence {
    String sequence

    // Constructor
    DNASequence(String sequence) {
        this.sequence = sequence.toUpperCase() // Ensure sequence is in uppercase for consistency
    }

    // Method to calculate melting temperature using the Wallace rule
    double getMeltingTemperature() {
        int g_count = sequence.count('G')
        int c_count = sequence.count('C')
        int a_count = sequence.count('A')
        int t_count = sequence.count('T')

        // Wallace rule calculation
        double tm = 4 * (g_count + c_count) + 2 * (a_count + t_count)
        return tm
    }

    String toString() {
        return "DNA[$sequence]"
    }
}
```

The `DNASequence` class is available in the execution context:

```nextflow
// main.nf
workflow {
    channel.of('ACGTTGCAATGCCGTA', 'GCGTACGGTACGTTAC')
        .map { seq -> new DNASequence(seq) }
        .view { dna ->
            "Found sequence '$dna' with melting temperaure ${dna.getMeltingTemperature()}°C"
        }
}
```

It prints:

```
Found sequence 'DNA[ACGTTGCAATGCCGTA]' with melting temperaure 48.0°C
Found sequence 'DNA[GCGTACGGTACGTTAC]' with melting temperaure 50.0°C
```

:::{note}
Package declarations in the `lib` directory are ignored. The package of a class is determined by the directory structure within the `lib` directory.

For example, if the above example were defined in `lib/utils/DNASequence.groovy`, the class would need to be referenced in pipeline scripts as `utils.DNASequence`.
:::

### Data

In general, input data should be provided by external sources using parameters which can be controlled by the user. This way, a pipeline can be easily reused to process different datasets which are appropriate for the pipeline.

Parameters can be declared with default values in the main script or in the configuration file:

```groovy
params.my_input = 'default input file'
params.my_output = 'default output path'
params.my_flag = false
// ...
```

When launching a pipeline, parameter values can be provided on the command line or in a params file (using the `-params-file` option). Options prefixed with a double dash (`--`) are interpreted as parameters:

```bash
nextflow run <your pipeline> --my_input /path/to/input/file --my_output /other/path --my_flag true
```

When a pipeline requires some small data that rarely changes, it may be easier to include the data in the pipeline repository. You can reference this data from the pipeline script in a portable manner (i.e. without relying on an absolute path) by using the `projectDir` implicit variable, which refers to the local copy of the pipeline repository.

The following example references the file `dataset/sequences.fa` in the pipeline repository:

```nextflow
sequences = file("$projectDir/dataset/sequences.fa")
sequences.splitFasta {
    println it
}
```

### Environment

The "environment" refers to any other aspects of the environment in which your pipeline is executed, such as environment variables and resource managers like SLURM.

Any environment variable that may be required by the tools in your pipeline can be defined under the `env` scope in your Nextflow configuration. For example:

```groovy
env {
  DELTA = 'hello'
  GAMMA = 'world'
}
```

Similarly, if you use an HPC scheduler like SLURM or a cloud batch service like AWS Batch to execute tasks in a distributed manner, you can use a configuration profile to define the settings for a given environment.

See {ref}`config-page` for more information about Nextflow configuration and {ref}`executor-page` for more information about executors.
