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

import com.google.gson.Gson;
import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import io.seqera.http.HxClient;
import org.eclipse.jgit.lib.*;
import org.eclipse.jgit.transport.Connection;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.UnsupportedEncodingException;
import java.net.URI;
import java.net.URLEncoder;
import java.net.http.HttpClient;
import java.net.http.HttpRequest;
import java.net.http.HttpResponse;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;
import java.util.stream.StreamSupport;

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
    private final Gson gson = new Gson();


    public SeqeraBaseConnection(TransportSeqera transport) {
        this.transport = transport;
        // Parse URI: seqera://bucket/path
        String dataLinkName = transport.getURI().getHost();
        this.repoPath = transport.getURI().getPath().substring(1); // Remove leading slash
        SeqeraGitCredentialsProvider credentials = (SeqeraGitCredentialsProvider) transport.getCredentialsProvider();
        this.endpoint = credentials.getEndpoint();
        this.workspaceId = credentials.getWorkspaceId();
        this.httpClient = createHttpClient(credentials);
        this.dataLink = findDataLink(dataLinkName);
        loadRefsMap();
        log.debug("Created Seqera Connection for dataLink={} path={}", dataLink.id, repoPath);
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
        byte[] content = downloadFile(path);
        if (content != null) {
            return new String(content);
        }
        return null;
    }

    /**
     * Find DatalinkId from the bucket name
     */
    protected DataLink findDataLink(String linkName) {
        try {
            String url = endpoint + "/data-links?status=AVAILABLE&search=" + linkName;
            if (workspaceId != null) {
                url += "&workspaceId=" + workspaceId;
            }

            HttpRequest request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .GET()
                .build();

            HttpResponse<String> response = httpClient.send(request, HttpResponse.BodyHandlers.ofString());

            if (response.statusCode() < 200 && response.statusCode() >= 300) {
                String message = response.body();
                if (message == null) message = "";
                message = message + " - HTTP " + response.statusCode();
                throw new RuntimeException("Failed to find data-link: " + linkName + message);
            }

            JsonObject responseJson = gson.fromJson(response.body(), JsonObject.class);
            JsonArray files = responseJson.getAsJsonArray("dataLinks");

            List<DataLink> dataLinks = java.util.stream.StreamSupport.stream(files.spliterator(), false)
                .map(JsonElement::getAsJsonObject)
                .filter(obj -> obj.get("name").getAsString().equals(linkName))
                .map(obj -> getDataLinkFromJSONObject(obj))
                .toList();

            if (dataLinks.isEmpty()) {
                log.debug("No datalink response: " + linkName);
                throw new RuntimeException("No Data-link found for " + linkName);
            }
            return dataLinks.getFirst();
        } catch (IOException e) {
            throw new RuntimeException("Exception finding data-link", e);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new RuntimeException("Download interrupted", e);
        }
    }

    private DataLink getDataLinkFromJSONObject(JsonObject obj) {
        String id = obj.get("id").getAsString();
        //Get first credential
        String credentialsId = obj.get("credentials").getAsJsonArray().get(0).getAsJsonObject().get("id").getAsString();

        return new DataLink(id, credentialsId);
    }

    /**
     * Download a file from the data-link
     */
    protected byte[] downloadFile(String filePath) throws IOException {
        try {
            Map<String, String> queryParams = new HashMap<>();
            if (workspaceId != null) {
                queryParams.put("workspaceId", workspaceId);
            }
            if (dataLink.credentialsId != null) {
                queryParams.put("credentialsId", dataLink.credentialsId);
            }

            String url = buildDataLinkUrl("/download/"+ URLEncoder.encode(filePath, StandardCharsets.UTF_8), queryParams);

            HttpRequest request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .GET()
                .build();

            HttpResponse<byte[]> response = httpClient.send(request, HttpResponse.BodyHandlers.ofByteArray());

            if (response.statusCode() == 404) {
                log.debug("File {} in data-link {} not found", filePath, dataLink.id);
                return null;
            }

            if (response.statusCode() >= 400) {
                throw new IOException("Failed to download file: " + filePath + " - status: " + response.statusCode());
            }

            return response.body();
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Download interrupted", e);
        }
    }

    /**
     * Download a file from the data-link directly to a local file
     * This is more memory-efficient for large files like bundles
     */
    protected void downloadFileToPath(String filePath, Path targetPath) throws IOException {
        try {
            Map<String, String> queryParams = new HashMap<>();
            if (workspaceId != null) {
                queryParams.put("workspaceId", workspaceId);
            }
            if (dataLink.credentialsId != null) {
                queryParams.put("credentialsId", dataLink.credentialsId);
            }

            String url = buildDataLinkUrl("/download/"+ URLEncoder.encode(filePath, StandardCharsets.UTF_8), queryParams);

            HttpRequest request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .GET()
                .build();

            HttpResponse<Path> response = httpClient.send(request, HttpResponse.BodyHandlers.ofFile(targetPath));

            if (response.statusCode() == 404) {
                log.debug("File {} in data-link {} not found", filePath, dataLink.id);
                throw new IOException("File not found: " + filePath);
            }

            if (response.statusCode() >= 400) {
                throw new IOException("Failed to download file: " + filePath + " - status: " + response.statusCode());
            }

            log.debug("Downloaded {} to {}", filePath, targetPath);
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Download interrupted", e);
        }
    }

    /**
     * Upload a file to the data-link
     */
    protected void uploadFile(String filePath, Path localFile) throws IOException {
        try {
            // Step 1: Get upload URL
            Map<String, String> queryParams = new HashMap<>();
            if (workspaceId != null) {
                queryParams.put("workspaceId", workspaceId);
            }
            if (dataLink.credentialsId != null) {
                queryParams.put("credentialsId", dataLink.credentialsId);
            }
            String url = buildDataLinkUrl("/upload", queryParams);

            JsonObject uploadRequest = new JsonObject();
            uploadRequest.addProperty("fileName", filePath);

            HttpRequest request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .header("Content-Type", "application/json")
                .POST(HttpRequest.BodyPublishers.ofString(gson.toJson(uploadRequest)))
                .build();

            HttpResponse<String> response = httpClient.sendAsString(request);

            if (response.statusCode() >= 400) {
                throw new IOException("Failed to get upload URL for: " + filePath + " - status: " + response.statusCode());
            }

            JsonObject responseJson = gson.fromJson(response.body(), JsonObject.class);
            String uploadUrl = responseJson.get("uploadUrl").getAsString();

            // Step 2: Upload file to the provided URL
            HttpRequest uploadFileRequest = HttpRequest.newBuilder()
                .uri(URI.create(uploadUrl))
                .header("Content-Type", "application/octet-stream")
                .PUT(HttpRequest.BodyPublishers.ofFile(localFile))
                .build();

            HttpResponse<Void> uploadResponse = httpClient.send(uploadFileRequest, HttpResponse.BodyHandlers.discarding());

            if (uploadResponse.statusCode() >= 400) {
                throw new IOException("Failed to upload file: " + filePath + " - status: " + uploadResponse.statusCode());
            }

            // Step 3: Finish upload (required for some providers like AWS)

            String finishUrl = buildDataLinkUrl("/upload/finish", Map.of());
            if (workspaceId != null) {
                finishUrl += "?workspaceId=" + workspaceId;
            }

            JsonObject finishRequest = new JsonObject();
            finishRequest.addProperty("fileName", filePath);

            HttpRequest finishReq = HttpRequest.newBuilder()
                .uri(URI.create(finishUrl))
                .header("Content-Type", "application/json")
                .POST(HttpRequest.BodyPublishers.ofString(gson.toJson(finishRequest)))
                .build();

            httpClient.send(finishReq, HttpResponse.BodyHandlers.discarding());

        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Upload interrupted", e);
        }
    }

    /**
     * List files in a path
     */
    protected List<String> listFiles(String path) throws IOException {
        Map<String, String> queryParams = new HashMap<>();
        if (workspaceId != null) {
            queryParams.put("workspaceId", workspaceId);
        }
        if (dataLink.credentialsId != null) {
            queryParams.put("credentialsId", dataLink.credentialsId);
        }
        String url = buildDataLinkUrl("/browse/" + URLEncoder.encode(path, StandardCharsets.UTF_8), queryParams);
        log.debug(" GET {}", url);
        HttpRequest request = HttpRequest.newBuilder()
            .uri(URI.create(url))
            .GET()
            .build();

        HttpResponse<String> response = httpClient.sendAsString(request);

        if (response.statusCode() == 404) {
            log.debug("No files found for {}", url);
            return List.of();
        }

        if (response.statusCode() >= 400) {
            log.debug("Error getting files {}", url);
            String message = response.body();
            if (message == null) message = "";
            message = message + " - HTTP " + response.statusCode();

            throw new IOException("Failed to list files in: " + path + message);
        }
        log.debug(response.body());
        JsonObject responseJson = gson.fromJson(response.body(), JsonObject.class);
        JsonArray files = responseJson.getAsJsonArray("objects");

        return StreamSupport.stream(files.spliterator(), false)
            .map(JsonElement::getAsJsonObject)
            .map(obj -> obj.get("name").getAsString())
            .collect(java.util.stream.Collectors.toList());
    }

    /**
     * Delete a file from the data-link
     */
    protected void deleteFile(String filePath) throws IOException {
        try {
            Map<String, String> queryParams = new HashMap<>();

            if (workspaceId != null) {
                queryParams.put("workspaceId", workspaceId);
            }
            String url = buildDataLinkUrl("/content", queryParams);

            JsonObject deleteRequest = new JsonObject();
            deleteRequest.addProperty("path", filePath);

            HttpRequest request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .header("Content-Type", "application/json")
                .method("DELETE", HttpRequest.BodyPublishers.ofString(gson.toJson(deleteRequest)))
                .build();

            HttpResponse<Void> response = httpClient.send(request, HttpResponse.BodyHandlers.discarding());

            if (response.statusCode() >= 400 && response.statusCode() != 404) {
                throw new IOException("Failed to delete file: " + filePath + " - status: " + response.statusCode());
            }
        } catch (InterruptedException e) {
            Thread.currentThread().interrupt();
            throw new IOException("Delete interrupted", e);
        }
    }

    private String buildDataLinkUrl(String path, Map<String, String> queryParams) {
        final StringBuilder url = new StringBuilder(endpoint + "/data-links/" + URLEncoder.encode(dataLink.id, StandardCharsets.UTF_8));
        if (path != null) {
            if (!path.startsWith("/")) url.append("/");
            url.append(path);
        }

        if (queryParams != null && !queryParams.isEmpty()) {
            String queryString = queryParams.entrySet().stream()
                .map(entry -> URLEncoder.encode(entry.getKey(), StandardCharsets.UTF_8)
                    + "="
                    + URLEncoder.encode(entry.getValue(), StandardCharsets.UTF_8))
                .collect(Collectors.joining("&"));
            url.append('?').append(queryString);
        }
        return url.toString();
    }

    private void loadRefsMap() {
        log.debug("Loading refs Maps");
        try {
            addRefs("heads");
        } catch (Exception e) {
            log.debug("No heads found for dataLink={} path={}: {}", dataLink.id, repoPath, e.getMessage());
        }
        try {
            addRefs("tags");
        } catch (Exception e) {
            log.debug("No tags found for dataLink={} path={}: {}", dataLink.id, repoPath, e.getMessage());
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
            log.debug("No default refs found for dataLink={} path={}: {}", dataLink.id, repoPath, e.getMessage());
        }
    }

    private void addRefs(String refType) throws IOException {
        String path = repoPath + "/refs/" + refType;
        List<String> branches = listFiles(path);

        if (branches == null || branches.isEmpty()) {
            log.debug("No {} refs found for dataLink={} path={}", refType, dataLink.id, repoPath);
            return;
        }

        for (String branch : branches) {
            String branchPath = path + "/" + branch;
            List<String> bundles = listFiles(branchPath);
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

    static class DataLink {
        private String id;
        private String credentialsId;

        private DataLink(String id, String credentialsId) {
            this.id = id;
            this.credentialsId = credentialsId;
        }

        public String getId() {
            return id;
        }

        public String getCredentialsId() {
            return credentialsId;
        }

    }
}
