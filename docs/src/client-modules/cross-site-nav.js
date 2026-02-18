(function () {
  if (typeof window === 'undefined') return;
  /**
   * Derive the product base from the first path segment at load time.
   * e.g. /nextflow/some/page → /nextflow/
   * Captured once so it doesn't drift as client-side nav happens.
   */
  const firstSegment = window.location.pathname.split('/').filter(Boolean)[0];
  const PRODUCT_BASE = firstSegment ? `/${firstSegment}/` : '/';

  function isExternalToProduct(rawUrl) {
    let url;
    try {
      url = new URL(rawUrl, window.location.href);
    } catch {
      return false;
    }

    // Ignore different origins — browser handles those natively
    if (url.origin !== window.location.origin) return false;

    // Force hard load if the path leaves the current product base
    return !url.pathname.startsWith(PRODUCT_BASE);
  }

  function hardNavigate(url) {
    window.location.assign(url instanceof URL ? url.href : url);
  }

  // ── 1. Intercept clicks before React Router sees them ───────────────────────
  document.addEventListener(
    'click',
    (e) => {
      const anchor = e.target.closest('a');
      if (!anchor) return;

      const href = anchor.getAttribute('href');
      if (!href || href.startsWith('#')) return;

      if (isExternalToProduct(anchor.href)) {
        e.preventDefault();
        e.stopImmediatePropagation();
        hardNavigate(anchor.href);
      }
    },
    true // capture phase — runs before React Router's synthetic event handling
  );

  // ── 2. Patch history.pushState — covers programmatic navigation ─────────────
  const _pushState = history.pushState.bind(history);
  history.pushState = function (state, title, url) {
    if (url && isExternalToProduct(String(url))) {
      hardNavigate(String(url));
      return;
    }
    return _pushState(state, title, url);
  };

  // ── 3. Patch history.replaceState — defensive, covers edge cases ────────────
  const _replaceState = history.replaceState.bind(history);
  history.replaceState = function (state, title, url) {
    if (url && isExternalToProduct(String(url))) {
      hardNavigate(String(url));
      return;
    }
    return _replaceState(state, title, url);
  };
})();
