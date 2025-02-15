(vscode-page)=

# VS Code integration

The [Nextflow VS Code extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) provides IDE support for Nextflow pipelines.

## Features

### Syntax highlighting

The VS Code extension highlights Nextflow scripts and config files for better readability.

### Diagnostics

The extension highlights source code in red for errors and yellow for warnings.

To view all diagnostics for the workspace, open the **Problems** tab. Here, you can search for diagnostics by diagnostic message, filename, and so on.

:::{note}
The language server parses scripts and config files according to the {ref}`Nextflow language specification <syntax-page>`, which is more strict than the Nextflow CLI. See {ref}`updating-syntax-page` for more information.
:::

### Hover hints

When you hover over certain source code elements, such as variable names and function calls, the extension provides a tooltip with related information, such as the definition and/or documentation for the element.

### Code navigation

The **Outline** section in the Explorer panel lists top-level definitions when you view a script. Include declarations in scripts and config files act as links, and ctrl-clicking them opens the corresponding script or config file.

To view the definition of a symbol (e.g., a workflow, process, function, or variable), right-click the symbol and select **Go to Definition**. Ctrl-click the symbol to view the definition. Ctrl-click the definition to show all references.

### Code completion

The extension suggests auto-completions for variable names, function names, config settings, and other symbols as you type. The extension also provides several snippets for common script declarations, such as processes and workflows.

### Formatting

The extension can format your scripts and config files based on a standard set of formatting rules. Rules can be customized using the **Nextflow > Formatting** [extension settings](#settings).

Use the **Format Document** command in the command palette to format the current file.

### Renaming symbols

The extension can rename all references of a symbol (e.g., a workflow, process, function, or variable) throughout the code.

To rename a symbol, right-click the symbol, select **Rename Symbol**, and enter a new name.

### Parameter schema

If a `nextflow_schema.json` file exists in the same directory as a script with an entry workflow, the extension uses the schema to provide validation, hover hints, and code completion for params in the entry workflow.

### DAG preview for workflows

The extension can generate a workflow DAG that includes the workflow inputs, outputs, and any processes or workflows that are called by the selected workflow. The workflow DAG is displayed in a new panel to the side.

To preview the DAG of a workflow, select the **Preview DAG** CodeLens above the workflow definition.

:::{note}
The **Preview DAG** CodeLens is only available when the script does not contain any errors.
:::

## Troubleshooting

In the event of an error, you can stop or restart the language server from the command palette. See [Commands](#commands) for the set of available commands.

Report issues at [nextflow-io/vscode-language-nextflow](https://github.com/nextflow-io/vscode-language-nextflow) or [nextflow-io/language-server](https://github.com/nextflow-io/language-server). When reporting, include a minimal code snippet that reproduces the issue and any error logs from the server. To view logs, open the **Output** tab and select **Nextflow Language Server** from the dropdown. Enable **Nextflow > Debug** in the [extension settings](#settings) to show additional log messages while debugging.

## Limitations

- The language server does not detect certain filesystem changes, such as changing the current Git branch. Restart the language server from the command palette to sync it with your workspace.

- The language server does not recognize configuration options from third-party plugins and will report "Unrecognized config option" warnings for them.

- The language server provides limited support for Groovy scripts in the `lib` directory. Errors in Groovy scripts are not reported as diagnostics, and changing a Groovy script does not automatically re-compile the Nextflow scripts that reference it. Edit the Nextflow script or close and re-open it to refresh the diagnostics.

## Commands

The following commands are available from the command palette:

- Restart language server
- Stop language server

(vscode-settings)=

## Settings

The following settings are available:

`nextflow.debug`
: Enable debug logging and debug information in hover hints.

`nextflow.files.exclude`
: Configure glob patterns for excluding folders from being searched for Nextflow scripts and configuration files.

`nextflow.formatting.harshilAlignment`
: Use the [Harshil Alignment™️](https://nf-co.re/docs/contributing/code_editors_and_styling/harshil_alignment) when formatting Nextflow scripts and config files.

`nextflow.java.home`
: Specifies the folder path to the JDK. Use this setting if the extension cannot find Java automatically.

`nextflow.paranoidWarnings`
: Enable additional warnings for future deprecations, potential problems, and other discouraged patterns.

## Language server

Most of the functionality of the VS Code extension is provided by the [Nextflow language server](https://github.com/nextflow-io/language-server), which implements the [Language Server Protocol (LSP)](https://microsoft.github.io/language-server-protocol/) for Nextflow scripts and config files.

The language server is distributed as a standalone Java application. It can be integrated with any editor that functions as an LSP client. Currently, only the VS Code integration is officially supported, but community contributions for other editors are welcome.
