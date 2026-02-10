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
          lastVersion: latest_version || undefined,
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
    }),
  });
}
