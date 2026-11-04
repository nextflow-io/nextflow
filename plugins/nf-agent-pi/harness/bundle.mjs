// Copyright 2013-2026, Seqera Labs
// SPDX-License-Identifier: Apache-2.0

// Bundles runner.mjs and its dependencies into a single file for the runner image, so the image
// ships no node_modules. Usage: node harness/bundle.mjs <outdir>

import { build } from "esbuild";
import { mkdirSync } from "node:fs";
import { join } from "node:path";

const outdir = process.argv[2];
if (!outdir) throw new Error("usage: node harness/bundle.mjs <outdir>");
mkdirSync(outdir, { recursive: true });

await build({
  // pi-ai loads its Bedrock provider through an import esbuild cannot follow, so register it
  // statically, as the Pi SDK does for its own Bun binary. The harness only reads models once a
  // start frame arrives, which is after this module body has run.
  stdin: {
    contents: [
      'import { bedrockProviderModule } from "@earendil-works/pi-ai/bedrock-provider";',
      'import { setBedrockProviderModule } from "@earendil-works/pi-ai/compat";',
      'import "./runner.mjs";',
      "setBedrockProviderModule(bedrockProviderModule);",
    ].join("\n"),
    resolveDir: import.meta.dirname,
    sourcefile: "entry.mjs",
  },
  outfile: join(outdir, "runner.mjs"),
  bundle: true,
  minify: true,
  platform: "node",
  format: "esm",
  target: "node24",
  // ESM output has no require or __dirname, which CommonJS dependencies still use
  banner: {
    js: [
      'import { createRequire as __nfCreateRequire } from "node:module";',
      'import { fileURLToPath as __nfFileURLToPath } from "node:url";',
      "const require = __nfCreateRequire(import.meta.url);",
      "const __filename = __nfFileURLToPath(import.meta.url);",
      'const __dirname = __nfFileURLToPath(new URL(".", import.meta.url));',
    ].join("\n"),
  },
  logLevel: "warning",
});
