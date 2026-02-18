# Versioning Documentation

This document explains how versioning is implemented for the Nextflow documentation.

## Overview

The Nextflow documentation uses Docusaurus versioning to maintain multiple versions of the documentation. This follows the same pattern as the platform-enterprise documentation in the unified docs repository.

## Files

- **`versions.json`**: Contains an array of version strings (e.g., `["24.10", "24.04"]`).
- **`latest_version.js`**: Exports the latest version from `versions.json` (first element in the array).
- **`docusaurus.config.js`**: Configuration file that includes versioning settings.

## Configuration

The versioning configuration in `docusaurus.config.js` includes:

- `includeCurrentVersion`: Set to `true` by default. Can be disabled with `INCLUDE_NEXT=false` environment variable for PR previews.
- `lastVersion`: Points to the latest version from `latest_version.js`, or `undefined` if no versions exist.

## Creating a New Version

To create a new version of the documentation, use the Docusaurus versioning command:

```bash
npm run docusaurus docs:version <version>
```

For example:

```bash
npm run docusaurus docs:version 24.10
```

This will:

1. Create a copy of the current docs in `versioned_docs/version-<version>/`
2. Create a copy of the current sidebar in `versioned_sidebars/version-<version>-sidebars.json`
3. Add the version to `versions.json`

## Version Naming

Nextflow uses a calendar-based versioning scheme:

- Stable releases: `XX.04` and `XX.10` (April and October)
- Edge releases: All other months (e.g., `XX.01`, `XX.02`, `XX.09`, etc.)

For documentation versions, use the format: `YY.MM` (e.g., `24.10`, `25.04`)

## Directory Structure

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

## Controlling Version Switcher Visibility

The version switcher will only appear when you have 2+ versions. This includes both:
- Versioned snapshots in `versions.json` (e.g., "24.10", "24.04")
- The "next" version (current docs folder) if `includeCurrentVersion: true`

**Examples:**
- `versions.json: ["24.10"]` + `includeCurrentVersion: true` = 2 versions → **Switcher shows**
- `versions.json: []` + `includeCurrentVersion: true` = 1 version → **No switcher**
- `versions.json: ["24.10", "24.04"]` + `includeCurrentVersion: true` = 3 versions → **Switcher shows**

## Environment Variables

- `INCLUDE_NEXT`: Set to `"false"` to exclude the current version (useful for PR previews)

## Version Selection

Users can select versions via a version switcher in the documentation sidebar.

**Important:** The version switcher only appears when there are **2 or more versions** available. If you only have one version, no switcher will be shown.

The version switcher is implemented as a custom theme component that:

- Appears in the sidebar when viewing `/nextflow/*` pages
- Shows the current version with a "(current)" label
- Lists all available versions from `versions.json`
- Preserves the current page path when switching versions (falls back to the main page if the current page doesn't exist in the selected version)
- Uses Tailwind CSS for responsive styling with light/dark mode support

### Custom Theme Components

Two theme components are swizzled (customized) from the Seqera theme:

1. **`src/theme/DocSidebar/Desktop/index.tsx`** - Parent component that controls when to show the version switcher
   - Modified to check for `/nextflow` paths instead of `/platform-enterprise`
   - Renders the VersionSwitcher component when on Nextflow docs pages

2. **`src/theme/DocSidebar/Desktop/Content/VersionSwitcher/index.tsx`** - The version switcher component itself
   - Uses `"default"` plugin ID instead of `"platform-enterprise"`
   - Checks for `/nextflow` paths
   - Uses Docusaurus hooks:
     - `useVersions("default")` - Gets all available versions
     - `useDocsVersion()` - Gets the current version
     - `useDocsPreferredVersion("default")` - Manages user's preferred version selection
