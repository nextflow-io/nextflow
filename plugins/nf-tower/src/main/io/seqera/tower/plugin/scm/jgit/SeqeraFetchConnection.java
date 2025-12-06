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
import org.eclipse.jgit.lib.Constants;
import org.eclipse.jgit.lib.ObjectId;
import org.eclipse.jgit.lib.ProgressMonitor;
import org.eclipse.jgit.lib.Ref;
import org.eclipse.jgit.transport.*;
import org.eclipse.jgit.util.FileUtils;

import java.io.IOException;
import java.io.OutputStream;
import java.net.URISyntaxException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.*;

import static io.seqera.tower.plugin.datalink.DataLinkUtils.downloadFile;
import static io.seqera.tower.plugin.datalink.DataLinkUtils.listFiles;

/**
 * Fetch Connection implementation for Seqera Platform data-links git-remote storage.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
public class SeqeraFetchConnection extends SeqeraBaseConnection implements FetchConnection {

    public SeqeraFetchConnection(TransportSeqera transport) {
        super(transport);
    }

    @Override
    public void fetch(ProgressMonitor monitor, Collection<Ref> want, Set<ObjectId> have) throws TransportException {
        Path tmpdir = null;
        try {
            tmpdir = Files.createTempDirectory("seqera-remote-git-");
            for (Ref r : want) {
                downloadBundle(r, have, tmpdir, monitor);
            }
        } catch (IOException e) {
            throw new TransportException(transport.getURI(), "Exception fetching branches", e);
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

    private void downloadBundle(Ref r, Set<ObjectId> have, Path tmpdir, ProgressMonitor monitor) throws IOException {
        log.debug("Fetching {} in {}", r.getName(), tmpdir);

        // List bundles in the ref directory
        String refPath = repoPath + "/" + r.getName();
        List<String> bundles = listFiles(httpClient, endpoint, dataLink, refPath, workspaceId);

        if (bundles == null || bundles.isEmpty()) {
            throw new TransportException(transport.getURI(), "No bundle for " + r.getName());
        }

        if (bundles.size() > 1) {
            throw new TransportException(transport.getURI(), "More than one bundle for " + r.getName());
        }

        String bundleName = bundles.get(0);
        String bundlePath = refPath + "/" + bundleName;

        Path localBundle = tmpdir.resolve(bundleName);
        Files.createDirectories(localBundle.getParent());
        log.trace("Downloading bundle {} for branch {} to {}", bundlePath, r.getName(), localBundle);

        // Download the bundle
        downloadFile( httpClient, endpoint, dataLink, bundlePath, localBundle, workspaceId);

        parseBundle(r, have, localBundle, monitor);
    }

    private void parseBundle(Ref r, Set<ObjectId> have, Path localBundle, ProgressMonitor monitor) throws TransportException {
        List<RefSpec> specs = new ArrayList<>();
        List<Ref> refs = new ArrayList<>();
        refs.add(r);
        specs.add(new RefSpec().setForceUpdate(true).setSourceDestination(Constants.R_REFS + '*', Constants.R_REFS + '*'));

        try (FetchConnection c = Transport.open(transport.getLocal(), new URIish(localBundle.toUri().toString())).openFetch(specs)) {
            c.fetch(monitor, refs, have);
        } catch (IOException | RuntimeException | URISyntaxException err) {
            close();
            throw new TransportException(transport.getURI(), err.getMessage(), err);
        }
    }

    @Override
    public void fetch(ProgressMonitor monitor, Collection<Ref> want, Set<ObjectId> have, OutputStream out) throws TransportException {
        fetch(monitor, want, have);
    }

    @Override
    public boolean didFetchIncludeTags() {
        return false;
    }

    @Override
    public boolean didFetchTestConnectivity() {
        return false;
    }

    @Override
    public void setPackLockMessage(String message) {
        // No pack lock message supported.
    }

    @Override
    public Collection<PackLock> getPackLocks() {
        return Collections.emptyList();
    }
}
