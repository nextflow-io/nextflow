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

import io.seqera.http.HxClient;
import io.seqera.tower.plugin.TowerHxClientFactory;
import io.seqera.tower.plugin.datalink.DataLink;
import nextflow.SysEnv;
import org.eclipse.jgit.lib.*;
import org.eclipse.jgit.transport.Connection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.net.http.HttpClient;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import static io.seqera.tower.plugin.datalink.DataLinkUtils.*;

/**
 * Base class for connections with Seqera Platform data-links git-remote compatibility
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 *
 */
public class SeqeraBaseConnection implements Connection {

    protected static final Logger log = LoggerFactory.getLogger(SeqeraBaseConnection.class);
    protected TransportSeqera transport;
    protected DataLink dataLink;
    protected String repoPath;
    protected HxClient httpClient;
    protected String endpoint;
    protected String workspaceId;
    private final Map<String, Ref> advertisedRefs = new HashMap<String, Ref>();


    public SeqeraBaseConnection(TransportSeqera transport) {
        this.transport = transport;
        // Parse URI: seqera://bucket/path
        String dataLinkName = transport.getURI().getHost();
        this.repoPath = transport.getURI().getPath().substring(1); // Remove leading slash
        SeqeraGitCredentialsProvider credentials = (SeqeraGitCredentialsProvider) transport.getCredentialsProvider();
        this.endpoint = credentials.getEndpoint();
        this.workspaceId = credentials.getWorkspaceId();
        this.httpClient = TowerHxClientFactory.httpClient(credentials.getAccessToken(), SysEnv.get("TOWER_REFRESH_TOKEN"), credentials.getEndpoint(), credentials.getRetryPolicy());
        this.dataLink = findDataLink(this.httpClient, this.endpoint, this.workspaceId, dataLinkName);
        loadRefsMap();
        log.debug("Created Seqera Connection for dataLink={} path={}", dataLink.getId(), repoPath);
    }

    protected HxClient createHttpClient(SeqeraGitCredentialsProvider credentials) {
        final HxClient.Builder builder = HxClient.newBuilder();

        // Set up authentication
        final String token = credentials.getAccessToken();
        // Count occurrences of '.' in token
        long dotCount = token.chars().filter(ch -> ch == '.').count();
        if (dotCount == 2) {
            builder.bearerToken(token);
        } else {
            try {
                final String plain = new String(java.util.Base64.getDecoder().decode(token));
                final int p = plain.indexOf('.');
                if (p != -1) {
                    builder.bearerToken(token);
                } else {
                    builder.basicAuth("@token:" + token);
                }
            } catch (Exception e) {
                builder.basicAuth("@token:" + token);
            }
        }

        return builder
            .followRedirects(HttpClient.Redirect.NORMAL)
            .version(HttpClient.Version.HTTP_1_1)
            .connectTimeout(java.time.Duration.ofSeconds(60))
            .build();
    }


    protected String getDefaultBranchRef() throws IOException {
        String path = repoPath + "/HEAD";
        byte[] content = getFileContent(httpClient, endpoint, dataLink, path, workspaceId);;
        if (content != null) {
            return new String(content);
        }
        return null;
    }

    private void loadRefsMap() {
        log.debug("Loading refs Maps");
        try {
            addRefs("heads");
        } catch (Exception e) {
            log.debug("No heads found for dataLink={} path={}: {}", dataLink.getId(), repoPath, e.getMessage());
        }
        try {
            addRefs("tags");
        } catch (Exception e) {
            log.debug("No tags found for dataLink={} path={}: {}", dataLink.getId(), repoPath, e.getMessage());
        }
        try {
            final String defaultBranch = getDefaultBranchRef();
            if (defaultBranch != null) {
                Ref target = advertisedRefs.get(defaultBranch);
                if (target != null) {
                    advertisedRefs.put(Constants.HEAD, new SymbolicRef(Constants.HEAD, target));
                }
            }
        } catch (Exception e) {
            log.debug("No default refs found for dataLink={} path={}: {}", dataLink.getId(), repoPath, e.getMessage());
        }
    }

    private void addRefs(String refType) throws IOException {
        String path = repoPath + "/refs/" + refType;
        List<String> branches = listFiles(httpClient, endpoint, dataLink, path, workspaceId);

        if (branches == null || branches.isEmpty()) {
            log.debug("No {} refs found for dataLink={} path={}", refType, dataLink.getId(), repoPath);
            return;
        }

        for (String branch : branches) {
            String branchPath = path + "/" + branch;
            List<String> bundles = listFiles(httpClient, endpoint, dataLink, branchPath, workspaceId);
            if (bundles != null && !bundles.isEmpty()) {
                // Get the first bundle (there should only be one per branch)
                String bundleName = bundles.get(0);
                addRef(refType, branch, bundleName);
            }
        }
    }

    private void addRef(String refType, String branchName, String bundleName) {
        String sha = bundleName.replace(".bundle", "");
        ObjectId objectId = ObjectId.fromString(sha);
        String refName = "refs/" + refType + "/" + branchName;
        if (refName.endsWith("/"))
            refName = refName.substring(0, refName.length() - 1);
        if ("heads".equals(refType)) {
            advertisedRefs.put(refName, new ObjectIdRef.PeeledNonTag(Ref.Storage.NETWORK, refName, objectId));
        } else if ("tags".equals(refType)) {
            advertisedRefs.put(refName, new ObjectIdRef.Unpeeled(Ref.Storage.NETWORK, refName, objectId));
        }
    }

    @Override
    public Map<String, Ref> getRefsMap() {
        return advertisedRefs;
    }

    @Override
    public Collection<Ref> getRefs() {
        return advertisedRefs.values();
    }

    @Override
    public Ref getRef(String name) {
        return advertisedRefs.get(name);
    }

    @Override
    public void close() {
        // Cleanup if needed
    }

    @Override
    public String getMessages() {
        return "";
    }

    @Override
    public String getPeerUserAgent() {
        return "";
    }

    static class BranchData {
        private String type;
        private String simpleName;
        private String refName;
        private ObjectId objectId;

        private BranchData(String type, String simpleName, ObjectId objectId) {
            this.type = type;
            this.simpleName = simpleName;
            this.refName = String.format("refs/%s/%s", type, simpleName);
            this.objectId = objectId;
        }

        public static BranchData fromKey(String key) {
            String[] parts = key.split("/");
            if (parts.length < 5) {
                throw new RuntimeException("Incorrect key format in Seqera git-remote repository. Key should include: repo-path/refs/<type>/<branch>/<hash>.bundle");
            }
            final String type = parts[parts.length - 3];
            final String branch = parts[parts.length - 2];
            final String sha = parts[parts.length - 1].replace(".bundle", "");
            return new BranchData(type, branch, ObjectId.fromString(sha));
        }

        public String getType() {
            return type;
        }

        public String getSimpleName() {
            return simpleName;
        }

        public String getRefName() {
            return refName;
        }

        public ObjectId getObjectId() {
            return objectId;
        }
    }
}
