import {
  createSeqeraConfig,
  getSeqeraThemeConfig,
  getSeqeraPresetOptions,
} from "@seqera/docusaurus-preset-seqera";

export default async function createConfigAsync() {
  return createSeqeraConfig({
    clientModules: [require.resolve('./src/client-modules/cross-site-nav.js')],
    themes: ["@docusaurus/theme-mermaid"],
    markdown: {
      mermaid: true,
    },
    plugins: [
      ['docusaurus-plugin-llms', {
        id: 'llms-nextflow',
        docsDir: 'docs',
        llmsTxtFilename: 'llms-nextflow.txt',
        llmsFullTxtFilename: 'llms-nextflow-full.txt',
        title: 'Nextflow',
        description: 'Documentation for Nextflow.',
        rootContent: 'This file contains links to Nextflow documentation following the llmstxt.org standard.',
        generateLLMsTxt: true,
        generateLLMsFullTxt: true,
        generateMarkdownFiles: true,
        includeBlog: false,
        excludeImports: true,
        removeDuplicateHeadings: true,
        ignoreFiles: ['**/tags', '**/tags/**'],
        processingBatchSize: 50,
      }],
      [
        "@docusaurus/plugin-content-docs",
        {
          id: "nextflow",
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
              host: 'uk4gflrza0d8yx5sp-1.a1.typesense.net',
              port: 443,
              protocol: 'https',
            },
          ],
          apiKey: 'KZsuSjc7jPqDm7pkl1kN8TkoHH9b3dwY',
          connectionTimeoutSeconds: 2,
        },
        typesenseSearchParameters: {
          query_by: 'content,hierarchy.lvl0,hierarchy.lvl1,hierarchy.lvl2,hierarchy.lvl3',
          group_by: 'url_without_anchor',
          query_by_weights: '1,1,4,3,2',
          drop_tokens_threshold: 0,
          group_limit: 1,
          num_typos: 1,
          prioritize_exact_match: true,
          filter_by: 'docusaurus_tag:!=[default,doc_tag_doc_list,blog_posts_list,blog_tags_posts,doc_tags_list,blog_tags_list,docs-platform-enterprise-25.2,docs-platform-enterprise-25.1,docs-platform-enterprise-24.2,docs-platform-enterprise-24.1,docs-platform-enterprise-23.4,docs-platform-enterprise-23.3]', // TODO Remove once the scraper is updated
        },
        contextualSearch: false,
        placeholder: 'Search Seqera docs...',
        productRoutes: [
          ['/platform-enterprise/', 'Platform Enterprise', 'platform-enterprise', null],
          ['/platform-cloud/', 'Platform Cloud', 'platform-cloud', null],
          ['/platform-cli/', 'Platform CLI', 'platform-cli', null],
          ['/platform-api/', 'Platform API', 'platform-api', null],
          ['/nextflow/', 'Nextflow', null, 'docs-nextflow-current'],
          ['/multiqc/', 'MultiQC', 'multiqc', null],
          ['/wave/', 'Wave', 'wave', null],
          ['/fusion/', 'Fusion', 'fusion', null],
          ['/changelog/', 'Changelog', null, null],
        ],
      },
      navbar: {
        items: [
          {
            label: 'Cloud',
            href: '/platform-cloud/',
            activeClassName: 'navbar__link--active',
          },
          {
            label: 'Enterprise',
            href: '/platform-enterprise/',
            activeClassName: 'navbar__link--active',
          },
          {
            label: 'Nextflow',
            href: '/nextflow/',
            activeClassName: 'navbar__link--active',
          },
          {
            label: 'MultiQC',
            href: '/multiqc/',
            activeClassName: 'navbar__link--active',
          },
          {
            label: 'Wave',
            href: '/wave/',
            activeClassName: 'navbar__link--active',
          },
          {
            label: 'Fusion',
            href: '/fusion/',
            activeClassName: 'navbar__link--active',
          },
        ],
      },
      seqera: {
        docs: {
          versionDropdown: {
            nextflow: {
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
