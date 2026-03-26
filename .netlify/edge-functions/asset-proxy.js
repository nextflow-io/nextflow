/**
 * Netlify Edge Function to proxy assets from multiple origins
 *
 * This function tries each origin in sequence until it finds the requested asset.
 * It handles both /assets/* and /img/* paths to ensure CSS, JavaScript, and images
 * are correctly served from the appropriate documentation site.
 *
 * The function checks origins in order and returns the first successful response.
 * If no origin has the asset, it returns a 404.
 */

export default async function handler(request, context) {
  const url = new URL(request.url);
  const assetPath = url.pathname; // e.g., /assets/css/styles.abc123.css

  // List of origins to try in order
  const origins = [
    'https://docs-migration.netlify.app',
    'https://seqera-docs-api.netlify.app',
    'https://seqera-docs.netlify.app',
    // Add more origins as needed
  ];

  // Try each origin in sequence
  for (const origin of origins) {
    try {
      const response = await fetch(`${origin}${assetPath}`);
      if (response.ok) {
        // Return the successful response directly
        return response;
      }
    } catch (error) {
      // Continue to next origin if fetch fails
      console.error(`Failed to fetch from ${origin}${assetPath}:`, error);
      continue;
    }
  }

  // If none of the origins worked, return 404
  return new Response('Asset not found', {
    status: 404,
    headers: {
      'Content-Type': 'text/plain',
    },
  });
}
