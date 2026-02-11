import {
  createSeqeraConfig,
  getSeqeraThemeConfig,
  getSeqeraPresetOptions,
} from "@seqera/docusaurus-preset-seqera";
import latest_version from "./latest_version.js";

export default async function createConfigAsync() {
  return createSeqeraConfig({
    themes: ["@docusaurus/theme-mermaid"],
    markdown: {
      mermaid: true,
    },
    plugins: [
      [
        "@docusaurus/plugin-content-docs",
        {
          routeBasePath: "/nextflow",
          path: "docs",
          sidebarPath: "./sidebars.js",
          showLastUpdateAuthor: false,
          showLastUpdateTime: false,
          // For PR Previews we want to see the latest doc-set with expected changes.
          includeCurrentVersion: process.env.INCLUDE_NEXT !== "false",
          // Set current as lastVersion to serve it at base path (/nextflow/)
          lastVersion: 'current',
          versions: {
            current: {
              label: 'Latest',
              banner: 'unreleased',
            },
            '24.10': {
              banner: 'none',
            },
          },
          rehypePlugins: [require("rehype-katex")],
          remarkPlugins: [
            (await import("remark-math")).default,
            (await import("remark-sectionize")).default,
          ],
        },
      ],
    ],
    presets: [
      [
        "@seqera/docusaurus-preset-seqera",
        await getSeqeraPresetOptions({
          docs: false,
          openapi: false,
          theme: {
            customCss: require.resolve("./src/custom.css"),
          },
        }),
      ],
    ],
    themeConfig: getSeqeraThemeConfig({
      prism: {
        additionalLanguages: ['nextflow', 'groovy', 'java', 'bash', 'yaml', 'json'],
      },
      navbar: {
        items: [
          {
            label: 'Cloud',
            href: '/platform-cloud/',
          },
          {
            label: 'Enterprise',
            href: '/platform-enterprise/',
          },
          {
            label: 'Nextflow',
            href: '/nextflow/',
          },
          {
            label: 'MultiQC',
            href: '/multiqc/',
          },
          {
            label: 'Wave',
            href: '/wave/',
          },
          {
            label: 'Fusion',
            href: '/fusion/',
          },
        ],
      },
      seqera: {
        docs: {
          versionDropdown: {
            default: {
              enabled: true,
              showCurrent: true,
              currentLabel: 'Latest',
            },
          },
        },
      },
    }),
  });
}
