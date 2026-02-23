# Nextflow documentation

This directory contains the official Nextflow documentation, built with [Docusaurus 3.9.2](https://docusaurus.io/) using the Seqera preset theme.

## Technology stack

- **Docusaurus 3.9.2**: Modern static site generator
- **MDX**: Markdown with JSX/React component support
- **Seqera Preset**: Custom Docusaurus theme and configuration
- **Node.js 20+**: Required runtime environment
- **Tailwind CSS**: Utility-first CSS framework
- **DaisyUI**: Component library

## Quick start

```bash
cd docs
npm install
npm start
```

This starts a local development server at `http://localhost:3000` with hot reload. Most changes are reflected live without restarting the server.

## Available commands

| Command | Description |
|---------|-------------|
| `npm start` | Start local development server with hot reload |
| `npm run build` | Build production static site to `build/` directory |
| `npm run serve` | Serve the built site locally for testing |
| `npm run clear` | Clear Docusaurus cache (use if build issues occur) |
| `npm run swizzle` | Eject Docusaurus components for customization |

## Directory structure

```
docs/
├── docs/                    # Documentation content (.mdx files)
│   ├── _static/            # Static assets for documentation pages
│   ├── developer/          # Developer/contributor documentation
│   ├── guides/             # How-to guides
│   ├── migrations/         # Version migration guides
│   ├── plugins/            # Plugin documentation
│   ├── reference/          # Language and runtime reference
│   ├── snippets/           # Code snippet examples
│   └── tutorials/          # Step-by-step tutorials
├── src/                    # Custom React components and styling
│   ├── components/         # Custom React components
│   └── custom.css          # Custom CSS styles
├── static/                 # Static assets served at root
├── build/                  # Generated static site (git-ignored)
├── docusaurus.config.js    # Main Docusaurus configuration
├── sidebars.js             # Sidebar navigation structure
└── package.json            # Node.js dependencies

../netlify.toml             # Netlify deployment configuration (at repo root)
```

## Writing documentation

### File format

- Use `.mdx` extension for all documentation files
- MDX supports standard Markdown plus JSX/React components
- Files are automatically processed and converted to HTML pages

### Navigation

- Add new pages to `sidebars.js` to appear in navigation
- Use relative paths for internal links: `[link text](./other-page.mdx)`
- Organize pages into categories using the sidebar structure

### Admonitions

Docusaurus provides built-in admonitions for callouts:

```markdown
:::note
This is a note
:::

:::tip
This is a helpful tip
:::

:::warning
This is a warning
:::

:::danger
This is a danger warning
:::
```

### Version tags

Use custom version tags for tracking feature changes:

```jsx
<AddedInVersion version="25.10.0" />
<ChangedInVersion version="25.10.0" />
<DeprecatedInVersion version="25.04.0" />
```

### Code blocks

Code blocks support syntax highlighting:

````markdown
```groovy
workflow {
    println "Hello, Nextflow!"
}
```
````

Supported languages include: `groovy`, `bash`, `python`, `javascript`, `java`, `yaml`, `json`, and many more.

### Images and assets

- Place images in `docs/_static/` or `static/`
- Reference from docs using relative paths:
  - From `docs/_static/`: `![alt text](./_static/image.png)`
  - From `static/`: `![alt text](/image.png)` (served at root)

## Versioning

This section explains how versioning is implemented.

### Files

- **`versions.json`**: Contains an array of version strings (e.g., `["24.10", "24.04"]`).
- **`latest_version.js`**: Exports the latest version from `versions.json` (first element in the array).
- **`docusaurus.config.js`**: Configuration file that includes versioning settings.

### Configuration

The versioning configuration in `docusaurus.config.js` includes:

- `includeCurrentVersion`: Set to `true` by default. Can be disabled with `INCLUDE_NEXT=false` environment variable for PR previews.
- `lastVersion`: Points to the latest version from `latest_version.js`, or `undefined` if no versions exist.

### Creating a new version

To create a new version of the documentation, use the Docusaurus versioning command:

```bash
npm run docusaurus docs:version <version>
```

This will:

1. Create a copy of the current docs in `versioned_docs/version-<version>/`
2. Create a copy of the current sidebar in `versioned_sidebars/version-<version>-sidebars.json`
3. Add the version to `versions.json`

### Directory structure

When versions are created, the structure will be:

```
docs/                           # Current/next version (latest documentation)
versioned_docs/
  version-24.10/               # October 2024 stable release docs
  version-24.04/               # April 2024 stable release docs
versioned_sidebars/
  version-24.10-sidebars.json  # Sidebar for 24.10
  version-24.04-sidebars.json  # Sidebar for 24.04
versions.json                   # ["24.10", "24.04"]
latest_version.js               # Exports "24.10"
```

### Controlling version switcher visibility

The version switcher will only appear when you have 2+ versions. This includes both:
- Versioned snapshots in `versions.json` (e.g., "24.10", "24.04")
- The "next" version (current docs folder) if `includeCurrentVersion: true`

#### Custom theme components

Two theme components are swizzled (customized) from the Seqera theme:

1. **`src/theme/DocSidebar/Desktop/index.tsx`** - Parent component that controls when to show the version switcher
   - Modified to check for `/nextflow` paths
   - Renders the VersionSwitcher component when on Nextflow docs pages

2. **`src/theme/DocSidebar/Desktop/Content/VersionSwitcher/index.tsx`** - The version switcher component itself
   - Uses `"default"` plugin ID instead of `"platform-enterprise"`
   - Checks for `/nextflow` paths
   - Uses Docusaurus hooks:
     - `useVersions("default")` - Gets all available versions
     - `useDocsVersion()` - Gets the current version
     - `useDocsPreferredVersion("default")` - Manages user's preferred version selection

## Contributing

We welcome documentation contributions! Please:

1. **Fork and create a branch** for your changes
2. **Test locally** using `npm run build`
3. **Follow existing patterns** for consistency
4. **Update sidebars.js** if adding new pages
5. **Check spelling and grammar**
6. **Submit a pull request** with a clear description

See the main [CONTRIBUTING.md](../CONTRIBUTING.md) for general contribution guidelines.

## Troubleshooting

### Build fails

1. Clear the cache: `npm run clear`
2. Reinstall dependencies: `rm -rf node_modules && npm install`
3. Check for syntax errors in MDX files
4. Verify all internal links are correct

### Hot reload not working

1. Restart the development server
2. Clear browser cache
3. Check for JavaScript errors in browser console

### Missing styling

1. Ensure custom CSS imports in `docusaurus.config.js`
2. Check Tailwind configuration
3. Clear Docusaurus cache

### Node version issues

Ensure you're using Node.js 20 or higher:

```bash
node --version
```

## License

Nextflow documentation is distributed under the [Creative Commons Attribution-ShareAlike 4.0 International (CC BY-SA 4.0)](https://creativecommons.org/licenses/by-sa/4.0/) license.