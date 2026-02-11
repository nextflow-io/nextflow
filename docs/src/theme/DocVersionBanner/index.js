import React from 'react';
import {useDocsVersion} from '@docusaurus/plugin-content-docs/client';

export default function DocVersionBanner() {
  const version = useDocsVersion();

  if (version.version !== 'current') {
    return null;
  }

  return (
    <div className="max-w-3xl mx-auto md:pr-4">
      <div className="theme-admonition theme-admonition-caution admonition_whTW alert alert--warning block shadow-none border border-gray-200">
        <div className="admonitionContent_pDMz">
          <p>
            The <strong>Latest</strong> Nextflow documentation is being migrated from nextflow.io. For the most up-to-date content, see{' '}
            <a href="https://nextflow.io/docs/latest/" target="_blank" rel="noopener noreferrer">
              <strong>https://nextflow.io/docs/</strong>
            </a>
            {' '}(latest).
          </p>
        </div>
      </div>
    </div>
  );
}
