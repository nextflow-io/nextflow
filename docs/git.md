(git-page)=

# Git

## Git configuration

The file `$HOME/.nextflow/scm` allows you to centralise the security credentials required to access private project repositories on Bitbucket, GitHub and GitLab source code management (SCM) platforms or to manage the configuration properties of private server installations (of the same platforms).

The configuration properties for each Git provider are defined inside the `providers` section. Properties for the same provider are grouped with a common name and delimited with curly brackets. For example:

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
: Git provider name, either: `github`, `gitlab` or `bitbucket`.

`providers.<provider>.server`
: *Required only for private SCM servers*
: SCM server name including the protocol prefix e.g. `https://github.com`.

`providers.<provider>.endpoint`
: *Required only for private SCM servers*
: SCM API `endpoint` URL e.g. `https://api.github.com` (default: the same as `providers.<provider>.server`).

## Git providers

### BitBucket

Create a `bitbucket` entry in the [SCM configuration file](#git-configuration) specifying your user name and app password, as shown below:

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

### BitBucket Server

[BitBucket Server](https://confluence.atlassian.com/bitbucketserver) is a self-hosted Git repository and management platform.

:::{note}
BitBucket Server uses a different API from the [BitBucket](https://bitbucket.org/) Cloud service. Make sure to use the right configuration whether you are using the cloud service or a self-hosted installation.
:::

To access your local BitBucket Server create an entry in the [SCM configuration file](#git-configuration) specifying as shown below:

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

### GitHub

Create a `github` entry in the [SCM configuration file](#git-configuration) specifying your user name and access token as shown below:

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

### GitLab

Create a `gitlab` entry in the [SCM configuration file](#git-configuration) specifying the user name, password and your API access token that can be found in your GitLab [account page](https://gitlab.com/profile/account) (sign in required). For example:

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

### Gitea

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

### Azure Repos

Nextflow has builtin support for [Azure Repos](https://azure.microsoft.com/en-us/services/devops/repos/), a Git source code management service hosted in the Azure cloud. To access your Azure Repos with Nextflow provide the repository credentials using the configuration snippet shown below:

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

### AWS CodeCommit

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
The `user` and `password` settings are optional. If omitted, the [AWS default credentials provider chain](https://docs.aws.amazon.com/sdk-for-java/v1/developer-guide/credentials.html) is used.
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

In order to use a private SCM installation you will need to set the server name and access credentials in your [SCM configuration file](#git-configuration) .

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

Or, alternatively, using the Git clone URL:

```bash
nextflow run http://gitlab.acme.org/foo/bar.git
```

:::{note}
You must also specify the server API endpoint URL if it differs from the server base URL. For example, for GitHub Enterprise V3, add `endpoint = 'https://git.your-domain.com/api/v3'`.
:::

:::{warning}
When accessing a private SCM installation over `https` from a server that uses a custom SSL certificate, you may need to import the certificate into your local Java keystore. See [Import the Certificate as a Trusted Certificate](https://docs.oracle.com/javase/tutorial/security/toolsign/rstep2.html) for more information.
:::
