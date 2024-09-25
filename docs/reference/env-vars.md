(config-env-vars)=

# Environment variables

The following environment variables control the configuration of the Nextflow runtime and the underlying Java virtual machine.

## Java settings

`JAVA_CMD`
: Defines the path location of the Java binary command used to launch Nextflow.

`JAVA_HOME`
: Defines the path location of the Java VM installation used to run Nextflow.

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
: For example, with `NXF_FILE_ROOT=/some/root/path`, the use of `file('foo')` will be resolved to the absolute path `/some/root/path/foo`. A remote root path can be specified using the usual protocol prefix, e.g. `NXF_FILE_ROOT=s3://my-bucket/data`. Files defined using an absolute path are not affected by this setting.

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

`NXF_PLUGINS_DEFAULT`
: Whether to use the default plugins when no plugins are specified in the Nextflow configuration (default: `true`).

`NXF_PLUGINS_DIR`
: The path where the plugin archives are loaded and stored (default: `$NXF_HOME/plugins`).

`NXF_PLUGINS_TEST_REPOSITORY`
: :::{versionadded} 23.04.0
  :::
: Defines a custom plugin registry or plugin release URL for testing plugins outside of the main registry. See {ref}`testing-plugins` for more information.

`NXF_PUBLISH_FAIL_ON_ERROR`
: :::{versionadded} 24.04.3
  :::
: Defines the default behavior of `publishDir.failOnError` setting. See {ref}`publishDir<process-publishdir>` directive for more information.

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

## Proxy settings

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

`NO_PROXY`
: Defines one or more host names that should not use the proxy server. Separate multiple names using a comma character.
