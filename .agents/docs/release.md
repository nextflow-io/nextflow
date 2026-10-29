# Release process

1. For each plugin changed since the last release, update its `CHANGELOG.md` and `VERSION` (semantic version by impact), and bump `nextflowVersion` in its `build.gradle` if it requires a newer Nextflow. Commit each plugin separately as `Bump <plugin>@<version>`, e.g. `Bump nf-amazon@2.0.0`.
   - On a stable branch (`STABLE-*`), the next semantic version has likely been published from master already, and the release would silently ship that build. Use a `-patchN` suffix on the last version released from the branch (e.g. `1.5.4-patch1`) and confirm it is unused at `https://registry.nextflow.io/api/v1/plugins/<name>`.
2. Update the root `VERSION`. Versions in April and October are stable (`25.10.0`); all other months are edge (`25.09.0-edge`).
3. Update the root `CHANGELOG.md` from `git log v<PREVIOUS VERSION>..`, leaving out `chore:` commits and docs-only changes (commits touching only `docs/`, `.docusaurus_site/`, or code comments).
4. Run `make check-agent-image` to confirm the nf-agent-pi runner image builds and its version matches its build context. On failure, bump the plugin `VERSION`, fix the Dockerfile, or both.
5. Run `make releaseInfo` to update the version number and generate checksums.
6. Stage the release files:
   ```
   git add \
     VERSION \
     CHANGELOG.md \
     nextflow \
     nextflow.md5 \
     nextflow.sha1 \
     nextflow.sha256 \
     modules/nextflow/src/main/resources/META-INF/plugins-info.txt \
     modules/nextflow/src/main/resources/META-INF/build-info.properties
   ```
7. Commit with the `[release]` and `[e2e prod]` tags and push to trigger the release workflow:
   ```
   git commit -s -m "[release] Nextflow version 25.09.0-edge [e2e prod]"
   git push origin master
   ```
