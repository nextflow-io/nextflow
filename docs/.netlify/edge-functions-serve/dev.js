import { boot } from "https://v2-17-1--edge.netlify.com/bootstrap/index-combined.ts";

const functions = {}; const metadata = { functions: {} };


      try {
        const { default: func } = await import("file:///Users/chris.hakkaart/workspace/architecture/nextflow/docs/netlify/edge-functions/asset-proxy.js");

        if (typeof func === "function") {
          functions["asset-proxy"] = func;
          metadata.functions["asset-proxy"] = {"url":"file:///Users/chris.hakkaart/workspace/architecture/nextflow/docs/netlify/edge-functions/asset-proxy.js"}
        } else {
          console.log("⬥ Failed to load Edge Function asset-proxy. The file does not seem to have a function as the default export.");
        }
      } catch (error) {
        console.log("⬥ Failed to run Edge Function asset-proxy:");
        console.error(error);
      }
      

boot(() => Promise.resolve(functions));