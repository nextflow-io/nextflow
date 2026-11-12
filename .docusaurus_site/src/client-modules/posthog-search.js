let debounceTimer;

function waitForPostHog(callback, maxAttempts = 50, interval = 100) {
  let attempts = 0;

  const checkPostHog = () => {
    if (window.posthog) {
      callback(window.posthog);
    } else if (attempts < maxAttempts) {
      attempts++;
      setTimeout(checkPostHog, interval);
    }
    // Silently fail if PostHog doesn't load - it's optional analytics
  };

  checkPostHog();
}

function setupTracking(ph) {
  document.addEventListener('input', (e) => {
    const input = e.target.closest('.DocSearch-Input, [class*="searchBox"] input');
    if (!input) return;

    clearTimeout(debounceTimer);
    debounceTimer = setTimeout(() => {
      const query = input.value;
      if (query && query.length > 2) {
        ph.capture('docs_search_performed', { query });
      }
    }, 500);
  });

  document.addEventListener('click', (e) => {
    const resultLink = e.target.closest('.DocSearch-Hit a, [class*="searchResult"] a');
    if (!resultLink) return;

    ph.capture('docs_search_result_clicked', {
      result_url: resultLink.href,
      result_title: resultLink.textContent?.trim() || '',
    });
  });
}

if (typeof window !== 'undefined') {
  waitForPostHog(setupTracking);
}