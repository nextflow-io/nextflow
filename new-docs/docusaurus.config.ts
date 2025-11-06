import type { Config } from "@docusaurus/types";
import { 
  createSeqeraConfig, 
  getSeqeraThemeConfig,
  getSeqeraPresetOptions 
} from "@seqera/docusaurus-preset-seqera";

export default async function createConfigAsync(): Promise<Config> {
  return createSeqeraConfig({
    plugins: [],
    presets: [
      [
        "@seqera/docusaurus-preset-seqera",
        getSeqeraPresetOptions({
        docs: {
          path: 'nextflow-docs',  // Path to YOUR docs folder
          routeBasePath: '/',  // or '/nextflow' if you want
          sidebarPath: './sidebars.ts',
        },
      }),
      ],
    ],

    themeConfig: getSeqeraThemeConfig({
 
    }),
  }) satisfies Config;
}