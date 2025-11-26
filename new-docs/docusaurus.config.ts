import type { Config } from "@docusaurus/types";
import {
  createSeqeraConfig,
  getSeqeraThemeConfig,
  getSeqeraPresetOptions,
} from "@seqeralabs/docusaurus-preset-seqera";

export default async function createConfigAsync(): Promise<Config> {
  return createSeqeraConfig({
    plugins: [],
    presets: [
      [
        "@seqeralabs/docusaurus-preset-seqera",
        await getSeqeraPresetOptions({  // Add 'await' here
          docs: {
            routeBasePath: "/nextflow",
            path: "nextflow-docs",
            sidebarPath: "./sidebars.ts",
            // editUrl: "https://github.com/nextflow/nextflow/tree/master/",
            showLastUpdateAuthor: false,
            showLastUpdateTime: false,
          },
          openapi: false,
          theme: {
            customCss: require.resolve("./src/css/custom.css"),
          },
        }),
      ],
    ],
    themeConfig: getSeqeraThemeConfig({}),
  }) satisfies Config;
}