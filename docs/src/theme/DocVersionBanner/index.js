import React from 'react';
import {useDocsVersion} from '@docusaurus/plugin-content-docs/client';

export default function DocVersionBanner() {
  const version = useDocsVersion();

  // Custom migration banner for current version
  if (version.version === 'current') {
    return (
      <div className="max-w-3xl mx-auto md:pr-4">
        <div className="theme-admonition theme-admonition-caution admonition_whTW alert alert--warning block shadow-none border border-gray-200">
          <div className="admonitionContent_pDMz">
            <p>
              Nextflow documentation is being migrated from{' '}
              <a href="https://nextflow.io/docs/" target="_blank" rel="noopener noreferrer">
                <strong>https://nextflow.io/docs/</strong>
              </a>
              {' '}. Both sites will remain active during the migration.
            </p>
          </div>
        </div>
      </div>
    );
  }

  // Default unreleased banner
  if (version.banner === 'unreleased') {
    return (
      <div className="max-w-3xl mx-auto md:pr-4">
        <div className="theme-admonition theme-admonition-caution admonition_whTW alert alert--warning block shadow-none border border-gray-200">
          <div className="admonitionContent_pDMz">
            <p>
              This is unreleased documentation for <strong>{version.label}</strong> version of Nextflow.
            </p>
          </div>
        </div>
      </div>
    );
  }

  // Default unmaintained banner
  if (version.banner === 'unmaintained') {
    return (
      <div className="max-w-3xl mx-auto md:pr-4">
        <div className="theme-admonition theme-admonition-caution admonition_whTW alert alert--warning block shadow-none border border-gray-200">
          <div className="admonitionContent_pDMz">
            <p>
              This is documentation for <strong>{version.label}</strong>, which is no longer actively maintained.
            </p>
          </div>
        </div>
      </div>
    );
  }

  return null;
}
