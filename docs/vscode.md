(vscode-page)=

# VS Code integration

The [Nextflow VS Code extension](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) provides language support for Nextflow pipelines in [VS Code](https://code.visualstudio.com/).

This page describes the [Nextflow language server](https://github.com/nextflow-io/language-server) used by the extension. See the extension README in [GitHub](https://github.com/nextflow-io/vscode-language-nextflow) or the [Visual Studio Marketplace](https://marketplace.visualstudio.com/items?itemName=nextflow.nextflow) for details on additional extension features.

## Features

### Syntax highlighting

The VS Code extension highlights Nextflow scripts and config files for better readability.

### Diagnostics

The extension highlights source code in red for errors and yellow for warnings.

To view all diagnostics for the workspace, open the **Problems** tab. Here, you can search for diagnostics by diagnostic message, filename, and so on.

:::{note}
The language server parses scripts and config files according to the {ref}`Nextflow language specification <syntax-page>`, which is more strict than the Nextflow CLI. See {ref}`strict-syntax-page` for more information.
:::

### Hover hints

When you hover over certain source code elements, such as variable names and function calls, the extension provides a tooltip with related information, such as the definition and/or documentation for the element.

If a [Javadoc](https://en.wikipedia.org/wiki/Javadoc) comment is defined above a workflow, process, or function, the extension will include the contents of the comment in hover hints. The following is an example Javadoc comment:

```nextflow
/**
 * Say hello to someone.
 *
 * @param name
 */
def sayHello(name) {
    println "Hello, ${name}!"
}
```

### Code navigation

The **Outline** section in the Explorer panel lists top-level definitions when you view a script. Include declarations in scripts and config files act as links, and ctrl-clicking them opens the corresponding script or config file.

To view the definition of a symbol (e.g., a workflow, process, function, or variable), right-click the symbol and select **Go to Definition**. Ctrl-click the symbol to view the definition. Ctrl-click the definition to show all references.

### Code completion

The extension suggests auto-completions for variable names, function names, config settings, and other symbols as you type. The extension also provides several snippets for common script declarations, such as processes and workflows.

### Formatting

The extension can format your scripts and config files based on a standard set of formatting rules. Rules can be customized using the **Nextflow > Formatting** extension settings.

Use the **Format Document** command in the command palette to format the current file.

### Renaming symbols

The extension can rename all references of a symbol (e.g., a workflow, process, function, or variable) throughout the code.

To rename a symbol, right-click the symbol, select **Rename Symbol**, and enter a new name.

### Parameter schema

If a `nextflow_schema.json` file exists in the same directory as a script with an entry workflow, the extension uses the schema to provide validation, hover hints, and code completion for params in the entry workflow.

### DAG preview for workflows

The extension can generate a workflow DAG that includes the workflow inputs, outputs, and any processes or workflows that are called by the selected workflow. The workflow DAG is displayed in a new panel to the side.

To preview the DAG of a workflow, select the **Preview DAG** CodeLens above the workflow definition.

## Troubleshooting

In the event of a language server error, you can use the **Nextflow: Restart language server** command in the command palette to restart the language server.

Report issues at [nextflow-io/vscode-language-nextflow](https://github.com/nextflow-io/vscode-language-nextflow) or [nextflow-io/language-server](https://github.com/nextflow-io/language-server). When reporting, include a minimal code snippet that reproduces the issue and any error logs from the server. To view logs, open the **Output** tab and select **Nextflow Language Server** from the dropdown. Enable **Nextflow > Debug** in the extension settings to show additional log messages while debugging.

## Limitations

- The language server does not detect certain filesystem changes, such as changing the current Git branch. Restart the language server from the command palette to sync it with your workspace.

- The language server does not recognize configuration options from third-party plugins and will report "Unrecognized config option" warnings for them.

- The language server provides limited support for Groovy scripts in the `lib` directory. Errors in Groovy scripts are not reported as diagnostics, and changing a Groovy script does not automatically re-compile the Nextflow scripts that reference it. Edit the Nextflow script or close and re-open it to refresh the diagnostics.

(vscode-language-server)=

## Language server

The Nextflow language server implements the [Language Server Protocol (LSP)](https://microsoft.github.io/language-server-protocol/) for Nextflow scripts and config files. It is distributed as a standalone Java application and can be integrated with any editor that functions as an LSP client.

Currently, only the VS Code integration is officially supported, but community contributions for other editors are welcome. Visit the [GitHub issues](https://github.com/nextflow-io/language-server/issues) page for the latest updates on community-led integrations.
