import React from 'react';
// Import the default MDX components from Docusaurus
import MDXComponents from '@theme-original/MDXComponents';
// Import the versioned admonition components from the Seqera theme
import {
  AddedInVersion,
  ChangedInVersion,
  DeprecatedInVersion,
} from '@seqera/docusaurus-theme-seqera/lib/theme/AdmonitionVersioned/AdmonitionVersioned';

export default {
  // Re-use the default mapping
  ...MDXComponents,
  // Add the versioned admonition components as global MDX components
  AddedInVersion,
  ChangedInVersion,
  DeprecatedInVersion,
};
