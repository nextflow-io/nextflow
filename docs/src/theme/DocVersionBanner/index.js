import React from 'react';
import {useDocsVersion} from '@docusaurus/plugin-content-docs/client';
import Admonition from '@theme/Admonition';

export default function DocVersionBanner() {
  const version = useDocsVersion();

  if (version.version !== 'current') {
    return null;
  }

  return (
    <div className="max-w-3xl mx-auto md:pr-4">
      <Admonition type="caution">
        <p>
          <strong>Nextflow documentation is currently being migrated from nextflow.io.</strong> For the latest content, see{' '}
          <a href="https://nextflow.io/docs/latest/" target="_blank" rel="noopener noreferrer">
            https://nextflow.io/docs/latest/
          </a>
        </p>
      </Admonition>
    </div>
  );
}
