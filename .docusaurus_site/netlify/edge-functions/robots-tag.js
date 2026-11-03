// Adds `X-Robots-Tag: noindex` for direct hits to docs-migration.netlify.app,
// while leaving responses unchanged when the request is proxied through
// docs.seqera.io (so the canonical URL is the only one indexed).
//
// Detection relies on the upstream proxy forwarding the original hostname in
// `X-Forwarded-Host`, which nginx, Cloudflare and most CDNs do by default.

export default async (request, context) => {
  const response = await context.next();

  const forwardedHost = request.headers.get("x-forwarded-host") || "";
  const isProxiedFromSeqera = forwardedHost.toLowerCase().includes("seqera.io");

  if (!isProxiedFromSeqera) {
    response.headers.set("X-Robots-Tag", "noindex");
  }

  return response;
};

export const config = {
  path: "/*",
};