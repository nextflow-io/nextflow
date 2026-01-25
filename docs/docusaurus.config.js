import {
  createSeqeraConfig,
  getSeqeraThemeConfig,
  getSeqeraPresetOptions,
} from "@seqera/docusaurus-preset-seqera";

export default async function createConfigAsync() {
  return createSeqeraConfig({
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
          remarkPlugins: [(await import("remark-math")).default],
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
    themeConfig: getSeqeraThemeConfig({}),
  });
}
