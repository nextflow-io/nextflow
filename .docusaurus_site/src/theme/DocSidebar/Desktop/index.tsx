import React, { useState } from 'react';
import clsx from 'clsx';
import {useThemeConfig} from '@docusaurus/theme-common';
import {useLocation} from '@docusaurus/router';
import Logo from '@theme/Logo';
import CollapseButton from '@theme/DocSidebar/Desktop/CollapseButton';
import Content from '@theme/DocSidebar/Desktop/Content';
import type {Props} from '@theme/DocSidebar/Desktop';
import VersionSwitcher from './Content/VersionSwitcher';
import styles from './styles.module.css';

function DocSidebarDesktop({path, sidebar, onCollapse, isHidden}: Props) {
  const {
    navbar: {hideOnScroll},
    docs: {
      sidebar: {hideable},
    },
  } = useThemeConfig();

  const location = useLocation();
  const [isOpen, setIsOpen] = useState(false);

  // Check if URL contains /nextflow
  const showVersionSwitcher = location.pathname.includes('/nextflow');

  return (
    <div
      className={clsx(
        styles.sidebar,
        hideOnScroll && styles.sidebarWithHideableNavbar,
        isHidden && styles.sidebarHidden,
        'h-full relative w-full'
      )}>
      {hideOnScroll && <Logo tabIndex={-1} className={styles.sidebarLogo} />}
      {showVersionSwitcher && (
        <VersionSwitcher isOpen={isOpen} setIsOpen={setIsOpen} />
      )}
      <Content path={path} sidebar={sidebar} />
      {hideable && <div className="absolute right-0 top-[3.75rem]"><CollapseButton onClick={onCollapse} /></div>}
    </div>
  );
}

export default React.memo(DocSidebarDesktop);
