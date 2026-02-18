import {
  createSeqeraConfig,
  getSeqeraThemeConfig,
  getSeqeraPresetOptions,
} from "@seqera/docusaurus-preset-seqera";
import latest_version from "./latest_version.js";

export default async function createConfigAsync() {
  return createSeqeraConfig({
    clientModules: [require.resolve('./src/client-modules/cross-site-nav.js')],
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
              banner: 'none',
            },
            '26.01.1-edge': {
              banner: 'none',
            },
            '24.10.4': {
              banner: 'none',
            },
          },
          rehypePlugins: [require("rehype-katex")],
          remarkPlugins: [
            (await import("remark-math")).default,
            (await import("remark-sectionize")).default,
            (await import("remark-code-snippets")).default,
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
      typesense: {
        typesenseCollectionName: 'seqera_docs',
        searchPagePath: '/search',
        typesenseServerConfig: {
          nodes: [
            {
              host: '9scwdgbn4v8r1lyfp.a1.typesense.net',
              port: 443,
              protocol: 'https',
            },
          ],
          apiKey: process.env.TYPESENSE_SEARCH_API_KEY,
          connectionTimeoutSeconds: 2,
        },
        typesenseSearchParameters: {
          query_by: 'content,hierarchy.lvl0,hierarchy.lvl1,hierarchy.lvl2,hierarchy.lvl3',
          group_by: 'url_without_anchor',
          group_limit: 1,
          num_typos: 1,
          prioritize_exact_match: true,
          filter_by: 'docusaurus_tag:!=[default,doc_tag_doc_list,blog_posts_list,blog_tags_posts,doc_tags_list,blog_tags_list]',
        },
        contextualSearch: false,
        placeholder: 'Search Seqera docs...',
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
