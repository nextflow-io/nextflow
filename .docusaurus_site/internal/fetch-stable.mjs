import git from "simple-git";
import fs from "fs-extra";
import path from "path";
import { fileURLToPath } from "url";
import { createRequire } from "module";

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);
const require = createRequire(import.meta.url);

// Fetch the latest stable tag from the remote (excludes edge)
async function getLatestStableTag(url) {
  const result = await git().listRemote(["--tags", "--refs", url]);
  const tags = result
    .split("\n")
    .map((line) => line.split("\t")[1])
    .filter((ref) => ref)
    .map((ref) => ref.replace("refs/tags/", ""))
    .filter((tag) => !tag.includes("edge"))
    .sort((a, b) => b.localeCompare(a, undefined, { numeric: true }));

  if (!tags.length) {
    throw new Error(`No stable tag found for ${url}`);
  }

  return tags[0];
}

const repositories = [
  {
    name: "nextflow",
    url: "https://github.com/nextflow-io/nextflow.git",
    path: path.join(__dirname, "..", "nextflow_stable"),
    getTag: (url) => getLatestStableTag(url),
  },
];

// Clone a repository at a specific ref
async function cloneRepo({ name, url, path: repoPath, ref }) {
  console.log(`${name}: Cloning into ${repoPath} at ${ref}`);
  try {
    await git().clone(url, repoPath, ["--branch", ref, "--depth", "1"]);
    console.log(`✅ ${name}: Cloned at ${ref}`);
  } catch (error) {
    console.error(`❌ Error cloning ${name}:`);
    console.error(error.message);
    process.exit(1);
  }
}

// Create versioned docs from a cloned repository
async function createVersionedDocs(repoPath, tag) {
  const srcDocs = path.join(repoPath, "docs", "docs");
  const srcSidebar = path.join(repoPath, "docs", "sidebars.js");
  const destDocs = path.join(__dirname, "..", "versioned_docs", `version-${tag}`);
  const destSidebar = path.join(__dirname, "..", "versioned_sidebars", `version-${tag}-sidebars.json`);
  const versionsFile = path.join(__dirname, "..", "versions.json");

  console.log(`Copying docs to ${destDocs}`);
  await fs.copy(srcDocs, destDocs, { overwrite: true });
  console.log(`✅ Docs copied`);

  console.log(`Converting sidebar to ${destSidebar}`);
  const sidebars = require(srcSidebar);
  await fs.ensureDir(path.dirname(destSidebar));
  await fs.writeJSON(destSidebar, sidebars, { spaces: 2 });
  console.log(`✅ Sidebar written`);

  await fs.writeJSON(versionsFile, [tag], { spaces: 2 });
  console.log(`✅ versions.json written: ${tag}`);
}

// Clone/update all repositories and create versioned docs
async function processRepositories() {
  const docsRoot = path.join(__dirname, "..");

  const toRemove = [
    ...repositories.map((repo) => repo.path),
    path.join(docsRoot, "versioned_docs"),
    path.join(docsRoot, "versioned_sidebars"),
    path.join(docsRoot, "versions.json"),
  ];

  for (const target of toRemove) {
    if (await fs.pathExists(target)) {
      console.log(`Removing ${target}`);
      await fs.remove(target);
    }
  }

  for (const repo of repositories) {
    try {
      const ref = repo.getTag ? await repo.getTag(repo.url) : repo.branch;
      await cloneRepo({ ...repo, ref });
      await createVersionedDocs(repo.path, ref);
    } catch (error) {
      console.error(`❌ Error processing ${repo.name}:`, error);
    }
  }
}

// Execute
processRepositories().catch(console.error);