(config-env-vars)=

# Environment variables

The following environment variables control the configuration of the Nextflow runtime and the underlying Java virtual machine.

## Java settings

`JAVA_CMD`
: Defines the path location of the Java binary command used to launch Nextflow.

`JAVA_HOME`
: Defines the path location of the Java VM installation used to run Nextflow.

(nxf-env-vars)=

## Nextflow settings

`NXF_ANSI_LOG`
: Enables/disables ANSI console output (default `true` when ANSI terminal is detected).

`NXF_ANSI_SUMMARY`
: Enables/disables ANSI completion summary: `true\|false` (default: print summary if execution last more than 1 minute).

`NXF_ASSETS`
: Defines the directory where downloaded pipeline repositories are stored (default: `$NXF_HOME/assets`)

`NXF_CACHE_DIR`
: :::{versionadded} 24.02.0-edge
  :::
: Defines the base cache directory when using the default cache store (default: `"$launchDir/.nextflow"`).

`NXF_CHARLIECLOUD_CACHEDIR`
: Directory where remote Charliecloud images are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_CLOUDCACHE_PATH`
: :::{versionadded} 23.07.0-edge
  :::
: Defines the base cache path when using the cloud cache store.

`NXF_CLOUD_DRIVER`
: Defines the default cloud driver to be used if not specified in the config file or as command line option, either `aws` or `google`.

`NXF_CONDA_CACHEDIR`
: Directory where Conda environments are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_CONDA_ENABLED`
: :::{versionadded} 22.08.0-edge
  :::
: Enable the use of Conda recipes defined by using the {ref}`process-conda` directive. (default: `false`).

`NXF_CONTAINER_ENTRYPOINT_OVERRIDE`
: :::{deprecated} 22.10.0
  :::
: When `true`, override the container entrypoint with `/bin/bash` (default: `false`).

`NXF_DATE_FORMAT`
: :::{versionadded} 25.07.0-edge
  :::
: Defines the format for date and time representations in notifications and reports. Supports custom formats (e.g., `yyyy-MM-dd HH:mm:ss`) or `iso` for ISO 8601 format with timezone (default: `dd-MMM-yyyy HH:mm:ss`).

`NXF_DEFAULT_DSL`
: :::{versionadded} 22.03.0-edge
  :::
: Defines the DSL version that should be used in not specified otherwise in the script of config file (default: `2`)

`NXF_DISABLE_CHECK_LATEST`
: :::{versionadded} 23.09.0-edge
  :::
: Nextflow automatically checks for a newer version of itself unless this option is enabled (default: `false`).

`NXF_DISABLE_JOBS_CANCELLATION`
: :::{versionadded} 21.12.0-edge
  :::
: Disables the cancellation of child jobs on workflow execution termination.

`NXF_DISABLE_PARAMS_TYPE_DETECTION`
: :::{versionadded} 23.07.0-edge
  :::
: Disables the automatic type detection of command line parameters.

`NXF_DISABLE_WAVE_SERVICE`
: :::{versionadded} 23.08.0-edge
  :::
: Disables the requirement for Wave service when enabling the Fusion file system.

`NXF_ENABLE_AWS_SES`
: :::{versionadded} 23.06.0-edge
  :::
: Enable to use of AWS SES native API for sending emails in place of legacy SMTP settings (default: `false`)

`NXF_ENABLE_FS_SYNC`
: :::{versionadded} 23.10.0
  :::
: When enabled the job script will execute Linux `sync` command on job completion. This may be useful to synchronize the job state over shared file systems (default: `false`)

`NXF_ENABLE_SECRETS`
: :::{versionadded} 21.09.0-edge
  :::
: Enable Nextflow secrets features (default: `true`)

`NXF_ENABLE_STRICT`
: :::{versionadded} 22.05.0-edge
  :::
: Enable Nextflow *strict* execution mode (default: `false`)

`NXF_ENABLE_VIRTUAL_THREADS`
: :::{versionadded} 23.05.0-edge
  :::
: Enable the use of virtual threads in the Nextflow runtime (default: `false`)

`NXF_EXECUTOR`
: Defines the default process executor, e.g. `sge`

`NXF_FILE_ROOT`
: :::{versionadded} 23.05.0-edge
  :::
: The file storage path against which relative file paths are resolved.
: For example, with `NXF_FILE_ROOT=/some/root/path`, the use of `file('hello')` will be resolved to the absolute path `/some/root/path/hello`. A remote root path can be specified using the usual protocol prefix, e.g. `NXF_FILE_ROOT=s3://my-bucket/data`. Files defined using an absolute path are not affected by this setting.

`NXF_HOME`
: Nextflow home directory (default: `$HOME/.nextflow`).

`NXF_JAVA_HOME`
: Defines the path location of the Java VM installation used to run Nextflow. This variable overrides the `JAVA_HOME` variable if defined.

`NXF_JVM_ARGS`
: :::{versionadded} 21.12.1-edge
  :::
: Allows the setting Java VM options. This is similar to `NXF_OPTS` however it's only applied the JVM running Nextflow and not to any java pre-launching commands.

`NXF_LOG_FILE`
: The filename of the Nextflow log (default: `.nextflow.log`)

`NXF_OFFLINE`
: When `true` prevents Nextflow from automatically downloading and updating remote project repositories (default: `false`).
: :::{versionchanged} 23.09.0-edge
  This option also disables the automatic version check (see `NXF_DISABLE_CHECK_LATEST`).
  :::
: :::{versionchanged} 23.11.0-edge
  This option also prevents plugins from being downloaded. Plugin versions must be specified in offline mode, or else Nextflow will fail.
  :::

`NXF_OPTS`
: Provides extra options for the Java and Nextflow runtime. It must be a blank separated list of `-Dkey[=value]` properties.

`NXF_ORG`
: Default `organization` prefix when looking for a hosted repository (default: `nextflow-io`).

`NXF_PARAMS_FILE`
: :::{versionadded} 20.10.0
  :::
: Defines the path location of the pipeline parameters file .

`NXF_PID_FILE`
: Name of the file where the process PID is saved when Nextflow is launched in background.

`NXF_PLUGINS_ALLOWED`
: :::{versionadded} 25.04.0
  :::
: Comma separated list of plugin IDs that can be used in a workflow executions e.g. `NXF_PLUGINS_ALLOWED=nf-amazon,nf-tower,nf-wave`. Use empty string to disallow all plugins.

`NXF_PLUGINS_DEFAULT`
: Whether to use the default plugins when no plugins are specified in the Nextflow configuration (default: `true`).

`NXF_PLUGINS_DIR`
: The path where the plugin archives are loaded and stored (default: `$NXF_HOME/plugins`).

`NXF_PLUGINS_REGISTRY_URL`
: :::{versionadded} 25.08.0-edge
  :::
: Specifies the URL of the plugin registry used to download and resolve plugins. This allows using custom or private plugin registries instead of the default public registry.

`NXF_PLUGINS_TEST_REPOSITORY`
: :::{versionadded} 23.04.0
  :::
: Defines a custom plugin registry or plugin release URL for testing plugins outside of the main registry.

`NXF_PUBLISH_FAIL_ON_ERROR`
: :::{versionadded} 24.04.3
  :::
: Defines the default behavior of `publishDir.failOnError` setting. See {ref}`publishDir<process-publishdir>` directive for more information.

`NXF_RETRY_POLICY_DELAY`
: :::{versionadded} 25.06.0-edge
  :::
: Delay used for HTTP retryable operations (default: `350ms`).

`NXF_RETRY_POLICY_JITTER`
: :::{versionadded} 25.06.0-edge
  :::
: Jitter value used for HTTP retryable operations (default: `0.25`).

`NXF_RETRY_POLICY_MAX_ATTEMPTS`
: :::{versionadded} 25.06.0-edge
  :::
: Max number of attempts used for HTTP retryable operations (default: `5`).

`NXF_RETRY_POLICY_MAX_DELAY`
: :::{versionadded} 25.06.0-edge
  :::
: Max delay used for HTTP retryable operations (default: `90s`).

`NXF_RETRY_POLICY_MULTIPLIER`
: :::{versionadded} 25.08.0-edge
  :::
: Delay multiplier used for HTTP retryable operations (default: `2.0`).

`NXF_SCM_FILE`
: :::{versionadded} 20.10.0
  :::
: Defines the path location of the SCM config file .

`NXF_SINGULARITY_CACHEDIR`
: Directory where remote Singularity images are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_SINGULARITY_LIBRARYDIR`
: :::{versionadded} 21.09.0-edge
  :::
: Directory where remote Singularity images are retrieved. It should be a directory accessible to all compute nodes.

`NXF_SPACK_CACHEDIR`
: Directory where Spack environments are stored. When using a computing cluster it must be a shared folder accessible from all compute nodes.

`NXF_SPACK_ENABLED`
: :::{versionadded} 23.02.0-edge
  :::
: Enable the use of Spack recipes defined by using the {ref}`process-spack` directive. (default: `false`).

`NXF_SYNTAX_PARSER`
: :::{versionadded} 25.02.0-edge
  :::
: Set to `'v2'` to use the {ref}`strict syntax <strict-syntax-page>` for Nextflow scripts and config files (default: `'v1'`).

`NXF_TEMP`
: Directory where temporary files are stored

`NXF_TRACE`
: Enable trace level logging for the specified packages. Equivalent to the `-trace` command-line option.

`NXF_VER`
: Defines which version of Nextflow to use.

`NXF_WORK`
: Directory where working files are stored (usually your *scratch* directory)

`NXF_WRAPPER_STAGE_FILE_THRESHOLD`
: :::{versionadded} 23.05.0-edge
  :::
: Defines the minimum size of the `.command.run` staging script for it to be written to a separate `.command.stage` file (default: `'1 MB'`).
: This setting is useful for executors that impose a size limit on job scripts.

## Seqera Platform settings

`TOWER_ACCESS_TOKEN`
: Access token for authenticating with Seqera Platform. Can also be configured via the `tower.accessToken` config option.

`TOWER_API_ENDPOINT`
: Defines the Seqera Platform API endpoint (default: `https://api.cloud.seqera.io`). Can also be configured via the `tower.endpoint` config option.

`TOWER_AUTH_DOMAIN`
: :::{versionadded} 25.10.0
  :::
: Specifies the Auth0 domain to use for authentication when connecting to a custom Platform endpoint. When set, this takes precedence over the built-in mappings for known Platform endpoints.

`TOWER_AUTH_ID`
: :::{versionadded} 25.10.0
  :::
: Specifies the Auth0 client ID to use for authentication when connecting to a custom Platform endpoint. Must be used in conjunction with `TOWER_AUTH_DOMAIN`.

`TOWER_REFRESH_TOKEN`
: Refresh token for maintaining authentication with Seqera Platform. Can also be configured via the `tower.refreshToken` config option.

`TOWER_WORKSPACE_ID`
: Workspace ID for the Seqera Platform workspace to use. Can also be configured via the `tower.workspaceId` config option.

## Other settings

`FTP_PROXY`
: :::{versionadded} 21.06.0-edge
  :::
: Defines the FTP proxy server. Proxy authentication is supported by providing the credentials in the proxy URL, e.g. `ftp://user:password@proxy-host.com:port`.

`HTTP_PROXY`
: Defines the HTTP proxy server.
: :::{versionadded} 21.06.0-edge
  Proxy authentication is supported by providing the credentials in the proxy URL, e.g. `http://user:password@proxy-host.com:port`.
  :::

`HTTPS_PROXY`
: Defines the HTTPS proxy server.
: :::{versionadded} 21.06.0-edge
  Proxy authentication is supported by providing the credentials in the proxy URL, e.g. `https://user:password@proxy-host.com:port`.
  :::

`NO_COLOR`
: Disables ANSI color codes in Nextflow log output. When this variable is set, Nextflow prints plain text logs following the [NO_COLOR standard](https://no-color.org/).
: If both `NO_COLOR` and `NXF_ANSI_LOG` are set, `NXF_ANSI_LOG` takes precedence.

`NO_PROXY`
: Defines one or more host names that should not use the proxy server. Separate multiple names using a comma character.

`TERMINAL_WIDTH`
: Forces the terminal width of ANSI-formatted log output. Overrides automatic terminal width detection and uses the specified width for line wrapping when set to an integer value.
