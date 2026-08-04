module.exports = {
    sidebar: [
        "index",
        {
            type: "category",
            label: "Getting started",
            collapsed: false,
            items: [
                "overview",
                "install",
                "developer-env",
                "your-first-script"
            ]
        },
        {
            type: "category",
            label: "Running pipelines",
            collapsed: true,
            items: [
                "cli",
                "config",
                "executor",
                "cache-and-resume",
                "reports"
            ]
        },
        {
            type: "category",
            label: "Developing pipelines",
            collapsed: true,
            items: [
                "script",
                "working-with-files",
                "process",
                "workflow",
                {
                    type: "category",
                    label: "Static typing",
                    collapsed: true,
                    link: { type: "doc", id: "static-typing" },
                    items: [
                        "process-typed",
                        "workflow-typed",
                        "typed-parameters"
                    ]
                },
                "notifications",
                "secrets",
                "sharing",
                "vscode"
            ]
        },
        {
            type: "category",
            label: "Modules",
            collapsed: true,
            link: { type: "doc", id: "modules/modules" },
            items: [
                "modules/using-modules",
                "modules/developing-modules",
                "modules/module-registry"
            ]
        },
        {
            type: "category",
            label: "Software dependencies",
            collapsed: true,
            items: [
                "git",
                "container",
                "conda",
                "spack",
                "wave"
            ]
        },
        {
            type: "category",
            label: "Compute and storage",
            collapsed: true,
            items: [
                "aws",
                "amazons3",
                "azure",
                "google",
                "kubernetes",
                "fusion"
            ]
        },
        {
            type: "category",
            label: "Plugins",
            collapsed: true,
            items: [
                "plugins/plugins",
                "plugins/using-plugins",
                "plugins/developing-plugins",
                "plugins/plugin-registry"
            ]
        },
        {
            type: "category",
            label: "Language reference",
            collapsed: true,
            items: [
                "reference/feature-flags",
                "reference/syntax",
                "reference/semantics",
                {
                    type: "category",
                    label: "Standard library",
                    link: { type: "doc", id: "reference/stdlib" },
                    items: [
                        {
                            type: "category",
                            label: "Namespaces",
                            link: { type: "doc", id: "reference/stdlib-namespaces" },
                            items: [
                                "reference/stdlib-namespaces/global",
                                "reference/stdlib-namespaces/channel",
                                "reference/stdlib-namespaces/log",
                                "reference/stdlib-namespaces/nextflow",
                                "reference/stdlib-namespaces/workflow"
                            ]
                        },
                        {
                            type: "category",
                            label: "Types",
                            link: { type: "doc", id: "reference/stdlib-types" },
                            items: [
                                "reference/stdlib-types/bag",
                                "reference/stdlib-types/boolean",
                                "reference/stdlib-types/channel",
                                "reference/stdlib-types/duration",
                                "reference/stdlib-types/float",
                                "reference/stdlib-types/integer",
                                "reference/stdlib-types/iterable",
                                "reference/stdlib-types/list",
                                "reference/stdlib-types/map",
                                "reference/stdlib-types/memory-unit",
                                "reference/stdlib-types/path",
                                "reference/stdlib-types/record",
                                "reference/stdlib-types/set",
                                "reference/stdlib-types/string",
                                "reference/stdlib-types/tuple",
                                "reference/stdlib-types/value",
                                "reference/stdlib-types/version-number"
                            ]
                        },
                        "reference/stdlib-groovy"
                    ]
                },
                "reference/process",
                "reference/channel",
                "reference/operator-typed",
                "reference/operator"
            ]
        },
        {
            type: "category",
            label: "Runtime reference",
            collapsed: true,
            items: [
                {
                    type: "category",
                    label: "CLI",
                    link: { type: "doc", id: "reference/cli" },
                    items: [
                        "reference/cli/auth",
                        "reference/cli/clean",
                        "reference/cli/clone",
                        "reference/cli/config",
                        "reference/cli/console",
                        "reference/cli/drop",
                        "reference/cli/fs",
                        "reference/cli/help",
                        "reference/cli/info",
                        "reference/cli/inspect",
                        "reference/cli/kuberun",
                        "reference/cli/launch",
                        "reference/cli/lineage",
                        "reference/cli/lint",
                        "reference/cli/list",
                        "reference/cli/log",
                        "reference/cli/logfile",
                        "reference/cli/module",
                        "reference/cli/plugin",
                        "reference/cli/pull",
                        "reference/cli/run",
                        "reference/cli/secrets",
                        "reference/cli/self-update",
                        "reference/cli/view"
                    ]
                },
                "reference/config",
                "reference/env-vars"
            ]
        },
        {
            type: "category",
            label: "Updates",
            collapsed: true,
            items: [
                "updating-nextflow",
                "strict-syntax",
                {
                    type: "category",
                    label: "Migration notes",
                    link: { type: "doc", id: "migrations/index" },
                    items: [
                        "migrations/26-04",
                        "migrations/25-10",
                        "migrations/25-04",
                        "migrations/24-10",
                        "migrations/24-04",
                        "migrations/dsl1"
                    ]
                }
            ]
        },
        {
            type: "category",
            label: "Contributing",
            collapsed: true,
            items: [
                "developer/index",
                "developer/diagram",
                "developer/config-scopes",
                {
                    type: "category",
                    label: "Packages",
                    link: { type: "doc", id: "developer/packages" },
                    items: [
                        "developer/nextflow",
                        "developer/nextflow.ast",
                        "developer/nextflow.cache",
                        "developer/nextflow.cli",
                        "developer/nextflow.cloud.aws",
                        "developer/nextflow.cloud.aws.nio",
                        "developer/nextflow.cloud.azure",
                        "developer/nextflow.cloud.google",
                        "developer/nextflow.config",
                        "developer/nextflow.container",
                        "developer/nextflow.dag",
                        "developer/nextflow.executor",
                        "developer/nextflow.extension",
                        "developer/nextflow.k8s",
                        "developer/nextflow.plugin",
                        "developer/nextflow.processor",
                        "developer/nextflow.scm",
                        "developer/nextflow.script",
                        "developer/nextflow.secret",
                        "developer/nextflow.trace"
                    ]
                }
            ]
        },
        {
            type: "category",
            label: "Tutorials",
            collapsed: true,
            items: [
                "tutorials/rnaseq-nf",
                "tutorials/data-lineage",
                "tutorials/workflow-outputs",
                "tutorials/static-types",
                "tutorials/static-types-operators",
                "tutorials/metrics",
                "tutorials/flux"
            ]
        },
        {
            type: "category",
            label: "Guides",
            collapsed: true,
            items: [
                "guides/aws-java-sdk-v2",
                "guides/gradle-plugin",
                "guides/migrate-plugin",
                "guides/updating-spot-retries"
            ]
        }
    ]
};
