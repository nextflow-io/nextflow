// Copyright 2013-2026, Seqera Labs
// SPDX-License-Identifier: Apache-2.0

// Entry point compiled into the runner image's harness executable (`bun build --compile`).
//
// pi-ai loads its Bedrock provider through a computed dynamic import that the bundler cannot
// follow, so register it statically, as the Pi SDK does for its own Bun binary. The harness only
// reads models once a start frame arrives, which is after this module body has run.

import { bedrockProviderModule } from "@earendil-works/pi-ai/bedrock-provider";
import { setBedrockProviderModule } from "@earendil-works/pi-ai/compat";
import "./runner.mjs";

setBedrockProviderModule(bedrockProviderModule);
