import {
  createSeqeraConfig,
  getSeqeraThemeConfig,
  getSeqeraPresetOptions,
} from "@seqera/docusaurus-preset-seqera";

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
