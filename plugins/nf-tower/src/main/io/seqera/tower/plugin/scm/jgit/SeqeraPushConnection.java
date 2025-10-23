/*
 * Copyright 2013-2025, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package io.seqera.tower.plugin.scm.jgit;

import org.eclipse.jgit.errors.TransportException;
import org.eclipse.jgit.lib.*;
import org.eclipse.jgit.revwalk.RevCommit;
import org.eclipse.jgit.revwalk.RevWalk;
import org.eclipse.jgit.transport.BundleWriter;
import org.eclipse.jgit.transport.PushConnection;
import org.eclipse.jgit.transport.RemoteRefUpdate;
import org.eclipse.jgit.util.FileUtils;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStream;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import java.util.Map;

/**
 * Push connection implementation for Seqera Platform data-links git-remote storage.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class SeqeraPushConnection extends SeqeraBaseConnection implements PushConnection {

    public SeqeraPushConnection(TransportSeqera transport) {
        super(transport);
    }

    @Override
    public void push(ProgressMonitor monitor, Map<String, RemoteRefUpdate> refUpdates) throws TransportException {
        Path tmpdir = null;
        try {
            tmpdir = Files.createTempDirectory("seqera-remote-git-");
            for (Map.Entry<String, RemoteRefUpdate> entry : refUpdates.entrySet()) {
                pushBranch(entry, tmpdir);
            }
        } catch (IOException e) {
            throw new TransportException(transport.getURI(), "Exception pushing branches", e);
        } finally {
            if (tmpdir != null) {
                try {
                    FileUtils.delete(tmpdir.toFile(), FileUtils.RECURSIVE);
                } catch (IOException e) {
                    throw new TransportException(transport.getURI(), "Exception cleaning up temporary files", e);
                }
            }
        }
    }

    private void pushBranch(Map.Entry<String, RemoteRefUpdate> entry, Path tmpdir) throws IOException {
        log.debug("Pushing {} reference", entry.getKey());
        final Ref ref = transport.getLocal().findRef(entry.getKey());
        if (ref == null || ref.getObjectId() == null) {
            throw new IllegalStateException("Branch " + entry.getKey() + " not found");
        }

        String oldBundlePath = checkExistingBundle(entry.getKey());
        if (oldBundlePath != null && isSameObjectId(oldBundlePath, ref.getObjectId())) {
            setUpdateStatus(entry.getValue(), RemoteRefUpdate.Status.UP_TO_DATE);
            return;
        }

        if (oldBundlePath != null && !isCommitInBranch(oldBundlePath, ref)) {
            setUpdateStatus(entry.getValue(), RemoteRefUpdate.Status.REJECTED_REMOTE_CHANGED);
            return;
        }

        log.trace("Generating bundle for branch {} in {}", entry.getKey(), tmpdir);
        Path bundleFile = bundle(ref, tmpdir);
        String bundlePath = repoPath + "/" + entry.getKey() + "/" + bundleFile.getFileName().toString();

        log.trace("Uploading bundle {} to data-link {}", bundleFile, dataLink.getId());
        uploadFile(bundlePath, bundleFile);

        if (oldBundlePath != null) {
            log.trace("Deleting old bundle {}", oldBundlePath);
            deleteFile(oldBundlePath);
        }

        setUpdateStatus(entry.getValue(), RemoteRefUpdate.Status.OK);

        if (getRef(Constants.HEAD) == null) {
            updateRemoteHead(entry.getKey());
        }
    }

    private void updateRemoteHead(String ref) {
        try {
            String headPath = repoPath + "/HEAD";
            // Try to download HEAD to check if it exists
            byte[] existingHead = downloadFile(headPath);
            if (existingHead == null) {
                log.debug("No remote default branch. Setting to {}.", ref);
                // Create a temporary file with the ref content
                Path tempFile = Files.createTempFile("head-", ".txt");
                try {
                    Files.write(tempFile, ref.getBytes(StandardCharsets.UTF_8));
                    uploadFile(headPath, tempFile);
                } finally {
                    Files.deleteIfExists(tempFile);
                }
            }
        } catch (IOException e) {
            log.warn("Failed to update remote HEAD", e);
        }
    }

    private boolean isSameObjectId(String bundlePath, ObjectId commitId) {
        BranchData branch = BranchData.fromKey(bundlePath);
        return branch.getObjectId().name().equals(commitId.name());
    }

    /**
     * Sets the status on a RemoteRefUpdate using reflection.
     * This is necessary because RemoteRefUpdate.status is package-private and JGit
     * doesn't provide a public API to set it. It also JAR signing verification which
     * disables the implementations of this class in the org.eclipse.jgit.transport package.
     * The Custom transport implementations like this Seqera transport need to update the
     * RemoteRefUpdate.status to inform callers of push results.
     */
    private void setUpdateStatus(RemoteRefUpdate update, RemoteRefUpdate.Status status) {
        try {
            Field statusField = RemoteRefUpdate.class.getDeclaredField("status");
            statusField.setAccessible(true);
            statusField.set(update, status);
        } catch (NoSuchFieldException e) {
            throw new RuntimeException("JGit API changed: RemoteRefUpdate.status field not found. " +
                "This may require updating the transport implementation.", e);
        } catch (IllegalAccessException e) {
            throw new RuntimeException("Unable to access RemoteRefUpdate.status field", e);
        } catch (Exception e) {
            throw new RuntimeException("Unexpected error setting status on RemoteRefUpdate", e);
        }
    }

    public boolean isCommitInBranch(String bundlePath, Ref branchRef) throws IOException {
        ObjectId commitId = BranchData.fromKey(bundlePath).getObjectId();
        try (RevWalk walk = new RevWalk(transport.getLocal())) {
            RevCommit branchTip = walk.parseCommit(branchRef.getObjectId());
            RevCommit targetCommit = walk.parseCommit(commitId);

            // Check if the commit is reachable from the branch tip
            return walk.isMergedInto(targetCommit, branchTip);
        }
    }

    private String checkExistingBundle(String refName) throws IOException {
        String refPath = repoPath + "/" + refName;
        List<String> bundles = listFiles(refPath);

        if (bundles == null || bundles.isEmpty()) {
            return null;
        }

        if (bundles.size() > 1) {
            throw new TransportException(transport.getURI(), "More than one bundle for " + refName);
        }

        return refPath + "/" + bundles.get(0);
    }

    @Override
    public void push(ProgressMonitor monitor, Map<String, RemoteRefUpdate> refUpdates, OutputStream out) throws TransportException {
        push(monitor, refUpdates);
    }

    private Path bundle(Ref ref, Path tmpdir) throws IOException {
        final BundleWriter writer = new BundleWriter(transport.getLocal());
        Path bundleFile = tmpdir.resolve(ref.getObjectId().name() + ".bundle");
        writer.include(ref);
        try (OutputStream out = new FileOutputStream(bundleFile.toFile())) {
            writer.writeBundle(NullProgressMonitor.INSTANCE, out);
        }
        return bundleFile;
    }

    @Override
    public void close() {
        // Cleanup if needed
    }
}
