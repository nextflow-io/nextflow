import React, { useState, useRef, useEffect, Dispatch, SetStateAction } from "react";
import {
  useVersions,
  useDocsVersion,
  useDocsPreferredVersion,
} from "@docusaurus/plugin-content-docs/client";
import { usePluginData } from "@docusaurus/useGlobalData";
import { useLocation } from "@docusaurus/router";
import Link from "@docusaurus/Link";

interface VersionSwitcherProps {
  isOpen: boolean;
  setIsOpen: Dispatch<SetStateAction<boolean>>;
}

const VersionSwitcher: React.FC<VersionSwitcherProps> = ({ isOpen, setIsOpen }) => {
  const dropdownRef = useRef<HTMLDivElement>(null);
  const location = useLocation();

  // check if the plugin exists first
  let pluginData;
  try {
    pluginData = usePluginData('docusaurus-plugin-content-docs', 'default');
  } catch (e) {
    return null;
  }

  const { savePreferredVersionName } = useDocsPreferredVersion("default");
  const versions = useVersions("default");
  const currentVersion = useDocsVersion();
  const [isNextflowPage, setNextflowPage] = useState(false);

  useEffect(() => {
    if (location.pathname.startsWith('/nextflow')) {
      setNextflowPage(true);
    } else {
      setNextflowPage(false);
    }
  }, [location.pathname]);

  useEffect(() => {
    const handleClickOutside = (event: MouseEvent) => {
      if (dropdownRef.current?.contains(event.target as Node)) return;
      setTimeout(() => setIsOpen(false), 100);
    };

    document.addEventListener("mousedown", handleClickOutside);
    return () => {
      document.removeEventListener("mousedown", handleClickOutside);
    };
  }, [setIsOpen]);

  const toggleDropdown = () => setIsOpen((prev) => !prev);

  function handleSelectVersion(version: string) {
    savePreferredVersionName(version);
  }

  // Get the display label for a version
  // Returns object with label and whether to show "v" prefix
  function getVersionLabel(version: any): { label: string; showV: boolean } {
    // Check if the label looks like a version number (contains dots or is all digits)
    const isVersionNumber = /^[\d.]+/.test(version.label);

    // Version numbers get "v" prefix, text labels like "Latest" don't
    return {
      label: version.label,
      showV: isVersionNumber
    };
  }

  if (typeof window === "undefined") return null;
  if (!isNextflowPage) return null;
  // Only show version switcher if there are multiple versions (2 or more)
  if (!versions || versions.length < 2) return null;

  const items = versions.filter(
    (version) => version.label !== currentVersion?.label,
  );

  let urlSuffix = "";

  if (
    currentVersion &&
    location.pathname.startsWith(`/nextflow/${currentVersion.label}`)
  ) {
    const currentVersionPrefix = `/nextflow/${currentVersion.label}`;
    urlSuffix = location.pathname.replace(currentVersionPrefix, "");
  }

  return (
    <div ref={dropdownRef} className="relative px-4">
      <div className={`bg-white dark:bg-gray-900 rounded-lg border border-gray-300 dark:border-gray-700 overflow-hidden ${isOpen ? 'rounded-b-none' : ''}`}>
        <button
          onClick={toggleDropdown}
          className="h-9 w-full flex items-center justify-between px-3 text-sm bg-transparent border-none cursor-pointer relative z-10 text-gray-900 dark:text-white transition-colors hover:bg-gray-100 dark:hover:bg-gray-800"
        >
          <span>
            {currentVersion ? (() => {
              const { label, showV } = getVersionLabel(currentVersion);
              return showV ? `v${label}` : label;
            })() : "Version"}
          </span>
          <svg
            className={`w-5 h-5 transition-transform duration-200 ${isOpen ? 'rotate-0' : '-rotate-90'}`}
            xmlns="http://www.w3.org/2000/svg"
            viewBox="4 4 16 16"
            fill="currentColor"
          >
            <path d="M11.8152 13.1989L10.0167 11.1432C9.80447 10.9013 9.97697 10.5214 10.2991 10.5214H13.8961C13.9682 10.5214 14.0388 10.5421 14.0994 10.5811C14.16 10.6201 14.2081 10.6758 14.2379 10.7414C14.2677 10.8071 14.2779 10.8799 14.2674 10.9512C14.2569 11.0226 14.226 11.0893 14.1785 11.1435L12.38 13.1985C12.3448 13.2388 12.3014 13.2711 12.2527 13.2932C12.204 13.3153 12.1511 13.3268 12.0976 13.3268C12.0441 13.3268 11.9912 13.3153 11.9425 13.2932C11.8938 13.2711 11.8504 13.2388 11.8152 13.1985V13.1989Z" />
          </svg>
        </button>

        {isOpen && (
          <div className="absolute left-4 right-4 top-full -mt-0.5 bg-white dark:bg-gray-900 border border-gray-300 dark:border-gray-700 border-t-gray-200 dark:border-t-gray-600 rounded-b-lg overflow-hidden z-50 transition-all duration-100">
            {items?.map((version) => (
              <div
                key={version.name}
                className="w-full"
                onClick={() => handleSelectVersion(version.name)}
              >
                <Link
                  to={`${version.path}${urlSuffix}`}
                  className="h-9 w-full flex items-center justify-between px-3 text-sm bg-transparent cursor-pointer text-gray-900 dark:text-white transition-colors hover:bg-gray-100 dark:hover:bg-gray-800 no-underline border-t border-gray-200 dark:border-gray-700 first:border-t-0"
                >
                  {(() => {
                    const { label, showV } = getVersionLabel(version);
                    return showV ? `v${label}` : label;
                  })()}
                </Link>
              </div>
            ))}
          </div>
        )}
      </div>
    </div>
  );
};

export default VersionSwitcher;
