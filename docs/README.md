# Nextflow documentation

This directory holds the Nextflow documentation content as `.mdx` files. The site is built with [Docusaurus](https://docusaurus.io/) and deployed to [docs.seqera.io/nextflow](https://docs.seqera.io/nextflow/) via Netlify.

## Editing docs

- Edit the `.mdx` files in this directory.
- Use relative links for internal references, e.g. `[link text](./other-page.mdx)`.
- Images and other assets live in `_static/`.
- Code snippets under `snippets/` are runnable and covered by the docs tests.

## Building and previewing locally

The Docusaurus project lives in [`.docusaurus_site/`](../.docusaurus_site) at the repository root (this directory is symlinked into it as `docs/`):

```bash
cd .docusaurus_site
npm install
npm start
```

See [`.docusaurus_site/README.md`](../.docusaurus_site/README.md) for the full build, configuration, and authoring guide.
