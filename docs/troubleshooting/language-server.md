(language-server-errors-page)=

# Language server errors

This page describes common language server errors and strategies to resolve them.

## Common errors

### Filesystem changes

The language server does not detect certain filesystem changes. For example, changing the current Git branch.

To resolve this issue, restart the language server from the command palette to sync it with your workspace. See [Stop and restart](#stop-and-restart) for more information.

### Third-party plugins

The language server does not recognize configuration options from third-party plugins and will report unrecognized config option warnings. There is currently no solution to suppress them.

### Groovy scripts

The language server provides limited support for Groovy scripts in the lib directory. Errors in Groovy scripts are not reported as diagnostics, and changing a Groovy script does not automatically re-compile the Nextflow scripts that reference it.

To resolve this issue, edit or close and re-open the Nextflow script to refresh the diagnostics.

## Stop and restart

In the event of an error, stop or restart the language server from the Command Palette. The following stop and restart commands are available:

- `Nextflow: Stop language server`
- `Nextflow: Restart language server`

See {ref}`vscode-commands` for a fill list of Nextflow VS Code extension commands.

## View logs

Error logs can be useful for troubleshooting errors.

To view logs in VS Code:

1. Open the **Output** tab in your console.
2. Select **Nextflow Language Server** from the dropdown.

To show additional log messages in VS Code:

1. Open the **Extensions** view in the left-hand menu.
2. Select the **Nextflow** extension.
3. Select the **Manage** icon.
3. Enable **Nextflow > Debug** in the extension settings.

## Report an issue

Report issues at [`nextflow-io/vscode-language-nextflow`](https://github.com/nextflow-io/vscode-language-nextflow) or [`nextflow-io/language-server`](https://github.com/nextflow-io/language-server). When reporting issues, include a minimal code snippet that reproduces the issue and any error logs from the server.
