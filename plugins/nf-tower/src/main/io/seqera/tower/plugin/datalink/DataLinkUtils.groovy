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

package io.seqera.tower.plugin.datalink

import groovy.json.JsonBuilder
import groovy.json.JsonSlurper
import groovy.transform.Canonical
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import io.seqera.http.HxClient

import java.net.http.HttpRequest
import java.net.http.HttpResponse
import java.nio.charset.StandardCharsets
import java.nio.file.Files
import java.nio.file.Path
import java.time.Duration

/**
 * Utility class for data-link file transfers (upload/download).
 * Implements cloud provider-specific upload strategies based on tower-cli patterns.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
@Slf4j
@CompileStatic
class DataLinkUtils {

    private static final JsonSlurper jsonSlurper = new JsonSlurper()
    private static final int CHUNK_SIZE = 250 * 1024 * 1024 // 250 MB chunks for multipart uploads
    private static final int DOWNLOAD_BUFFER_SIZE = 8192 // 8 KB buffer for downloads

    /**
     * Get file content as byte array from a data-link
     *
     * @param httpClient HTTP client to use
     * @param endpoint Platform API endpoint
     * @param dataLink Data-link associated to the file
     * @param filePath Remote file path
     * @param workspaceId Workspace ID (optional)
     * @return File content as byte array, or null if file not found
     * @throws IOException if download fails
     */
    static byte[] getFileContent(HxClient httpClient, String endpoint, DataLink dataLink,
                                 String filePath, String workspaceId) throws IOException {
        try {
            final credentialsId = dataLink.getCredentialsId()
            final dataLinkId = dataLink.getId()

            def queryParams = [:] as Map<String, String>
            if (workspaceId) {
                queryParams.workspaceId = workspaceId
            }
            if (credentialsId) {
                queryParams.credentialsId = credentialsId
            }

            def url = buildDataLinkUrl(endpoint, dataLinkId,
                "/download/${URLEncoder.encode(filePath, StandardCharsets.UTF_8)}", queryParams)

            def request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .GET()
                .build()

            def response = httpClient.send(request, HttpResponse.BodyHandlers.ofByteArray())

            if (response.statusCode() == 404) {
                log.debug("File {} in data-link {} not found", filePath, dataLinkId)
                return null
            }

            if (response.statusCode() >= 400) {
                throw new IOException("Failed to download file: ${filePath} - status: ${response.statusCode()}")
            }

            return response.body()

        } catch (InterruptedException e) {
            Thread.currentThread().interrupt()
            throw new IOException("Download interrupted", e)
        }
    }

    /**
     * Download a file from a data-link using streaming to avoid memory issues
     *
     * @param httpClient HTTP client to use
     * @param endpoint Platform API endpoint
     * @param dataLink Data-link associated to the file
     * @param filePath Remote file path
     * @param targetPath Local target path
     * @param workspaceId Workspace ID (optional)
     * @param credentialsId Credentials ID (optional)
     * @throws IOException if download fails
     */
    static void downloadFile(HxClient httpClient, String endpoint, DataLink dataLink,
                            String filePath, Path targetPath,
                            String workspaceId) throws IOException {
        try {
            final credentialsId = dataLink.getCredentialsId()
            final dataLinkId = dataLink.getId()
            // Build URL for getting download URL
            def queryParams = [:] as Map<String, String>
            if (workspaceId) {
                queryParams.workspaceId = workspaceId
            }
            if (credentialsId) {
                queryParams.credentialsId = credentialsId
            }
            queryParams.filePath = filePath
            queryParams.preview = 'false'

            def url = buildDataLinkUrl(endpoint, dataLinkId, "/generate-download-url", queryParams)

            def request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .timeout(Duration.ofSeconds(30))
                .GET()
                .build()

            def response = httpClient.sendAsString(request)
            log.debug("Upload response: {}", response.body())
            if (response.statusCode() == 404) {
                log.debug("File {} in data-link {} not found", filePath, dataLinkId)
                throw new IOException("File not found: ${filePath}")
            }

            if (response.statusCode() >= 400) {
                throw new IOException("Failed to get download URL for file: ${filePath} - status: ${response.statusCode()}")
            }

            // Parse response to get download URL
            def responseJson = jsonSlurper.parseText(response.body()) as Map
            def downloadUrl = responseJson.url as String

            // Download file using streaming
            log.debug("Downloading {} from {}", filePath, downloadUrl)
            downloadFileFromUrl(downloadUrl, targetPath)

        } catch (InterruptedException e) {
            Thread.currentThread().interrupt()
            throw new IOException("Download interrupted", e)
        }
    }

    /**
     * Download a file from a direct URL using streaming
     */
    private static void downloadFileFromUrl(String url, Path targetPath) throws IOException, InterruptedException {
        HxClient httpClient = HxClient.newHxClient()
        def request = HttpRequest.newBuilder()
            .uri(URI.create(url))
            .timeout(Duration.ofSeconds(30))
            .GET()
            .build()

        def response = httpClient.send(request, HttpResponse.BodyHandlers.ofInputStream())

        if (response.statusCode() != 200) {
            throw new IOException("Failed to download file: HTTP ${response.statusCode()}")
        }

        // Stream the download to avoid loading entire file in memory
        response.body().withCloseable { input ->
            Files.newOutputStream(targetPath).withCloseable { output ->
                def buffer = new byte[DOWNLOAD_BUFFER_SIZE]
                int bytesRead
                while ((bytesRead = input.read(buffer)) != -1) {
                    output.write(buffer, 0, bytesRead)
                }
            }
        }

        log.debug("Downloaded to {}", targetPath)
    }

    /**
     * Upload a file to a data-link with provider-specific strategies
     *
     * @param httpClient HTTP client to use
     * @param endpoint Platform API endpoint
     * @param dataLinkId Data-link ID
     * @param filePath Remote file path
     * @param localFile Local file to upload
     * @param provider Cloud provider type
     * @param workspaceId Workspace ID (optional)
     * @param credentialsId Credentials ID (optional)
     * @throws IOException if upload fails
     */
    static void uploadFile(HxClient httpClient, String endpoint, DataLink dataLink,
                          String filePath, Path localFile,
                          String workspaceId) throws IOException {
        try {
            if (!workspaceId)
                 throw new IOException("Workspace must be specified for uploading files to data links")
            final credentialsId = dataLink.getCredentialsId()
            final provider = dataLink.getProvider()
            final dataLinkId = dataLink.getId()
            // Step 1: Get upload URLs
            final queryParams = [:] as Map<String, String>
            queryParams.workspaceId = workspaceId
            if (credentialsId) {
                queryParams.credentialsId = credentialsId
            }

            def url = buildDataLinkUrl(endpoint, dataLinkId, "/upload", queryParams)
            // Detect MIME type
            def mimeType = Files.probeContentType(localFile) ?: 'application/octet-stream'

            def uploadRequest = [
                fileName: filePath,
                contentLength: localFile.toFile().length(),
                contentType: mimeType
            ]

            def request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .header('Content-Type', 'application/json')
                .POST(HttpRequest.BodyPublishers.ofString(new JsonBuilder(uploadRequest).toString()))
                .build()

            def response = httpClient.sendAsString(request)

            if (response.statusCode() >= 400) {
                def message = response.body() ?: ''
                message = "${message} - HTTP ${response.statusCode()}"
                throw new IOException("Failed to get upload URL for '${filePath}': ${message}")
            }

            log.debug("Upload response: {}", response.body())
            def responseJson = jsonSlurper.parseText(response.body()) as Map
            def uploadId = responseJson.uploadId as String
            def uploadUrls = responseJson.uploadUrls as List<String>

            // Step 2: Upload file using provider-specific strategy
            def etags = [] as List<String>
            def withError = false

            try {
                if (provider?.equalsIgnoreCase('AWS') || provider?.equalsIgnoreCase('SEQERACOMPUTE')) {
                    etags = uploadFileAws( localFile.toFile(), uploadUrls)
                } else if (provider?.equalsIgnoreCase('GOOGLE')) {
                    uploadFileGoogle( localFile.toFile(), uploadUrls)
                } else if (provider?.equalsIgnoreCase('AZURE')) {
                    uploadFileAzure( localFile.toFile(), uploadUrls)
                } else {
                    throw new IOException("Unsupported provider: ${provider}")
                }
            } catch (Exception e) {
                withError = true
                log.error("Upload failed: {}", e.message, e)
                throw new IOException("Failed to upload file: ${e.message}", e)
            } finally {
                // Step 3: Finalize upload
                finalizeUpload(httpClient, endpoint, dataLinkId, filePath, uploadId, etags, withError, workspaceId, credentialsId)
            }

        } catch (InterruptedException e) {
            Thread.currentThread().interrupt()
            throw new IOException("Upload interrupted", e)
        }
    }

    /**
     * Upload an entire folder (directory) to a data-link recursively
     *
     * @param httpClient HTTP client to use
     * @param endpoint Platform API endpoint
     * @param dataLink Data-link to upload to
     * @param remotePath Remote folder path in data-link
     * @param localFolder Local folder to upload
     * @param workspaceId Workspace ID
     * @throws IOException if upload fails
     */
    static void uploadFolder(HxClient httpClient, String endpoint, DataLink dataLink,
                            String remotePath, Path localFolder,
                            String workspaceId) throws IOException {
        if (!Files.isDirectory(localFolder)) {
            throw new IOException("Local path is not a directory: ${localFolder}")
        }

        log.debug("Uploading folder {} to {}", localFolder, remotePath)

        // Walk through the directory tree
        Files.walk(localFolder).each { Path localFile ->
            if (Files.isRegularFile(localFile)) {
                // Calculate relative path
                def relativePath = localFolder.relativize(localFile).toString()
                // Normalize path separators to forward slashes for remote path
                relativePath = relativePath.replace('\\', '/')
                // Combine with remote path
                def remoteFilePath = (remotePath ? "${remotePath}/${relativePath}" : relativePath) as String

                log.debug("Uploading file {} to {}", localFile, remoteFilePath)
                uploadFile(httpClient, endpoint, dataLink, remoteFilePath, localFile, workspaceId)
            }
        }

        log.debug("Folder upload completed: {}", remotePath)
    }

    /**
     * Download an entire folder (directory) from a data-link recursively
     *
     * @param httpClient HTTP client to use
     * @param endpoint Platform API endpoint
     * @param dataLink Data-link to download from
     * @param remotePath Remote folder path in data-link
     * @param localFolder Local folder to download to
     * @param workspaceId Workspace ID (optional)
     * @throws IOException if download fails
     */
    static void downloadFolder(HxClient httpClient, String endpoint, DataLink dataLink,
                              String remotePath, Path localFolder,
                              String workspaceId) throws IOException {
        log.debug("Downloading folder {} to {}", remotePath, localFolder)

        // Create local folder if it doesn't exist
        if (!Files.exists(localFolder)) {
            Files.createDirectories(localFolder)
        }

        // List all files in the remote folder
        def items = listFiles(httpClient, endpoint, dataLink, remotePath, workspaceId)

        for (String itemName : items) {
            // Calculate full remote path
            def remoteItemPath = (remotePath ? "${remotePath}/${itemName}" : itemName) as String
            def localItemPath = localFolder.resolve(itemName.replaceAll('/$', ''))

            // Check if item is a folder (ends with /)
            if (itemName.endsWith('/')) {
                // It's a folder, recurse
                log.debug("Downloading subfolder {}", remoteItemPath)
                def remoteFolderPath = remoteItemPath.replaceAll('/$', '') as String
                downloadFolder(httpClient, endpoint, dataLink, remoteFolderPath, localItemPath, workspaceId)
            } else {
                // It's a file, download it
                log.debug("Downloading file {} to {}", remoteItemPath, localItemPath)
                downloadFile(httpClient, endpoint, dataLink, remoteItemPath, localItemPath, workspaceId)
            }
        }

        log.debug("Folder download completed: {}", remotePath)
    }

    /**
     * Upload file using AWS S3 multipart upload protocol
     */
    private static List<String> uploadFileAws( File file, List<String> uploadUrls) throws IOException, InterruptedException {
        def etags = [] as List<String>
        int index = 0

        HxClient httpClient = HxClient.newHxClient()
        for (String url : uploadUrls) {
            log.debug("PUT chunk $index in $url")
            def chunk = getChunk(file, index)

            def request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .header('Content-Type', 'application/octet-stream')
                .PUT(HttpRequest.BodyPublishers.ofByteArray(chunk))
                .build()

            def response = httpClient.sendAsString(request)
            log.debug("Response: $response")
            if (response.statusCode() != 200) {
                log.error("Failed to upload chunk ${index}: HTTP ${response.statusCode()}, Message: ${response.body()}")
                throw new IOException("Failed to upload chunk ${index}: HTTP ${response.statusCode()}, Message: ${response.body()}")
            }

            // Extract ETag from response
            def etag = response.headers().firstValue('ETag').orElse(null)
            if (etag) {
                etags << etag
            } else {
                throw new IOException("Failed to get ETag from upload response - possible CORS issue")
            }
            log.debug("Chunk $index uploaded.")
            index++
        }

        return etags
    }

    /**
     * Upload file using Google Cloud Storage resumable upload protocol
     */
    private static void uploadFileGoogle( File file, List<String> uploadUrls) throws IOException, InterruptedException {
        if (!uploadUrls) {
            throw new IOException('No upload URLs provided')
        }

        def url = uploadUrls[0] // Google uses single resumable upload URL
        def fileSize = file.length()
        int index = 0
        long bytesUploaded = 0
        HxClient httpClient = HxClient.newHxClient()
        while (bytesUploaded < fileSize) {
            def chunk = getChunk(file, index)
            def start = bytesUploaded
            def end = Math.min(start + chunk.length, fileSize) - 1

            def request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .header('Content-Type', 'application/octet-stream')
                .header('Content-Range', "bytes ${start}-${end}/${fileSize}")
                .PUT(HttpRequest.BodyPublishers.ofByteArray(chunk))
                .build()

            def response = httpClient.sendAsString(request)

            if (response.statusCode() == 308) {
                // Resume incomplete - continue with next chunk
                def range = response.headers().firstValue('Range').orElse('')
                log.debug("Upload progress: {}", range)
            } else if (response.statusCode() == 200 || response.statusCode() == 201) {
                // Upload complete
                break
            } else {
                throw new IOException("Failed to upload chunk ${index}: HTTP ${response.statusCode()}")
            }

            bytesUploaded += chunk.length
            index++
        }
    }

    /**
     * Upload file using Azure Block Blob protocol
     */
    private static void uploadFileAzure(File file, List<String> uploadUrls) throws IOException, InterruptedException {
        def blockIds = [] as List<String>
        int index = 0
        HxClient httpClient = HxClient.newHxClient()
        for (String url : uploadUrls) {
            def chunk = getChunk(file, index)

            // Extract block ID from URL (typically in comp=block&blockid=XXX)
            def blockId = extractBlockId(url, index)
            blockIds << blockId

            def request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .header('Content-Type', 'application/octet-stream')
                .PUT(HttpRequest.BodyPublishers.ofByteArray(chunk))
                .build()

            def response = httpClient.sendAsString(request)

            if (response.statusCode() != 201 && response.statusCode() != 200) {
                throw new IOException("Failed to upload block ${index}: HTTP ${response.statusCode()}")
            }

            index++
        }
    }

    /**
     * Extract block ID from Azure upload URL
     */
    private static String extractBlockId(String url, int index) {
        // Try to extract from URL parameter
        def blockIdIdx = url.indexOf('blockid=')
        if (blockIdIdx != -1) {
            def endIdx = url.indexOf('&', blockIdIdx)
            if (endIdx == -1) endIdx = url.length()
            return url.substring(blockIdIdx + 8, endIdx)
        }
        // Fallback: generate block ID
        return String.format('block-%05d', index)
    }

    /**
     * Get a chunk of the file for multipart upload
     */
    private static byte[] getChunk(File file, int index) throws IOException {
        new RandomAccessFile(file, 'r').withCloseable { raf ->
            def start = (long) index * CHUNK_SIZE
            def end = Math.min(start + CHUNK_SIZE, file.length())
            def length = (int) (end - start)

            def buffer = new byte[length]
            raf.seek(start)
            raf.readFully(buffer)

            return buffer
        }
    }

    /**
     * Finalize the upload by calling the finish endpoint
     */
    private static void finalizeUpload(HxClient httpClient, String endpoint, String dataLinkId,
                                       String filePath, String uploadId, List<String> etags,
                                       boolean withError, String workspaceId, String credentialsId)
        throws IOException, InterruptedException {

        def queryParams = [:] as Map<String, String>
        if (workspaceId) {
            queryParams.workspaceId = workspaceId
        }
        if (credentialsId) {
            queryParams.credentialsId = credentialsId
        }

        def url = buildDataLinkUrl(endpoint, dataLinkId, "/upload/finish", queryParams)

        def finishRequest = [
            fileName: filePath,
            uploadId: uploadId,
            withError: withError
        ] as Map<String, Object>

        // Add ETags if present (for AWS)
        if (etags) {
            def tags = etags.withIndex().collect { String etag, int i ->
                [
                    eTag: etag,
                    partNumber: i + 1
                ] as Map<String, Object>
            }
            finishRequest.tags = tags
        }

        def request = HttpRequest.newBuilder()
            .uri(URI.create(url))
            .header('Content-Type', 'application/json')
            .POST(HttpRequest.BodyPublishers.ofString(new JsonBuilder(finishRequest).toString()))
            .build()

        def response = httpClient.sendAsString(request)

        if (response.statusCode() >= 400) {
            log.warn("Failed to finalize upload: HTTP {}, body: {}", response.statusCode(), response.body())
        } else {
            log.debug("Upload finalized successfully")
        }
    }

    /**
     * Build data-link URL with path and query parameters
     */
    private static String buildDataLinkUrl(String endpoint, String dataLinkId, String path, Map<String, String> queryParams) {
        def url = new StringBuilder("${endpoint}/data-links/${URLEncoder.encode(dataLinkId, StandardCharsets.UTF_8)}")
        if (path) {
            if (!path.startsWith('/')) url.append('/')
            url.append(path)
        }

        if (queryParams) {
            def queryString = queryParams.collect { key, value ->
                "${URLEncoder.encode(key, StandardCharsets.UTF_8)}=${URLEncoder.encode(value, StandardCharsets.UTF_8)}"
            }.join('&')
            url.append('?').append(queryString)
        }
        return url.toString()
    }

    /**
     * Find DatalinkId from the bucket name
     */
    static DataLink findDataLink(HxClient httpClient, String endpoint, String workspaceId, String linkName) {
        try {

            String url = endpoint + "/data-links?status=AVAILABLE&search=" + linkName
            if( workspaceId != null ) {
                url += "&workspaceId=" + workspaceId
            }

            HttpRequest request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .GET()
                .build()

            HttpResponse<String> response = httpClient.send(request, HttpResponse.BodyHandlers.ofString())

            if( response.statusCode() < 200 && response.statusCode() >= 300 ) {
                String message = response.body()
                if( message == null ) message = ""
                message = message + " - HTTP " + response.statusCode()
                throw new RuntimeException("Failed to find data-link: " + linkName + message)
            }

            def responseJson = jsonSlurper.parseText(response.body()) as Map
            def dataLinksArray = responseJson.dataLinks as List<Map>

            List<DataLink> dataLinks = dataLinksArray
                .findAll { it.name == linkName }
                .collect { getDataLinkFromJSON(linkName, it) }

            if( dataLinks.isEmpty() ) {
                log.debug("No datalink response: " + linkName)
                throw new RuntimeException("No Data-link found for " + linkName)
            }
            return dataLinks.first()
        } catch( IOException e ) {
            throw new RuntimeException("Exception finding data-link", e)
        } catch( InterruptedException e ) {
            Thread.currentThread().interrupt()
            throw new RuntimeException("Download interrupted", e)
        }
    }

    private static DataLink getDataLinkFromJSON(String linkName, Map obj) {
        String id = obj.id as String
        //Get first credential
        def credentials = obj.credentials as List<Map>
        String credentialsId = credentials[0].id as String
        //Get provider type (AWS, GOOGLE, AZURE, etc.)
        String provider = obj.provider ? obj.provider as String : "AWS"

        return new DataLink(linkName, id, credentialsId, provider)
    }

    /**
     * List files in a path within a data-link
     *
     * @param httpClient HTTP client to use
     * @param endpoint Platform API endpoint
     * @param dataLink Data-link to browse
     * @param path Path within the data-link to list
     * @param workspaceId Workspace ID (optional)
     * @return List of file/directory names in the specified path
     * @throws IOException if listing fails
     */
    static List<String> listFiles(HxClient httpClient, String endpoint, DataLink dataLink,
                                  String path, String workspaceId) throws IOException {
        try {
            final queryParams = [:] as Map<String, String>
            if (workspaceId) {
                queryParams.workspaceId = workspaceId
            }
            if (dataLink.getCredentialsId()) {
                queryParams.credentialsId = dataLink.getCredentialsId()
            }
            path = path ?: ''
            def url = buildDataLinkUrl(endpoint, dataLink.getId(),
                "/browse/${URLEncoder.encode(path, StandardCharsets.UTF_8)}", queryParams)
            log.debug("GET {}", url)

            def request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .GET()
                .build()

            def response = httpClient.sendAsString(request)

            if (response.statusCode() == 404) {
                log.debug("No files found for {}", url)
                return []
            }

            if (response.statusCode() >= 400) {
                log.debug("Error getting files {}", url)
                def message = response.body() ?: ''
                message = "${message} - HTTP ${response.statusCode()}"
                throw new IOException("Failed to list files in '${path}': ${message}")
            }
            
            def responseJson = jsonSlurper.parseText(response.body()) as Map
            def files = responseJson.objects as List<Map>

            return files.collect { it.name as String }

        } catch (InterruptedException e) {
            Thread.currentThread().interrupt()
            throw new IOException("List files interrupted", e)
        }
    }

    /**
     * Get the details of a path within a data-link
     *
     * @param httpClient HTTP client to use
     * @param endpoint Platform API endpoint
     * @param dataLink Data-link to browse
     * @param path Path within the data-link to get the details
     * @param workspaceId Workspace ID (optional)
     * @return Detail of the file/directory in the specified path. Null if file is not found
     * @throws IOException if listing fails
     */
    static DataLinkItem getFileDetails(HxClient httpClient, String endpoint, DataLink dataLink,
                                  String pathStr, String workspaceId) throws IOException {
        try {
            if(!pathStr)
                return new DataLinkItem(DataLinkItemType.FOLDER, dataLink.name, 0, null)

            def path = Path.of(pathStr)
            def fileName = path.getFileName()?.toString()
            def parent = path.getParent()?.toString() ?: ''
            final queryParams = [:] as Map<String, String>
            if (workspaceId) {
                queryParams.workspaceId = workspaceId
            }
            if (dataLink.getCredentialsId()) {
                queryParams.credentialsId = dataLink.getCredentialsId()
            }
            queryParams.search = fileName

            def url = buildDataLinkUrl(endpoint, dataLink.getId(),
                "/browse/${URLEncoder.encode(parent, StandardCharsets.UTF_8)}", queryParams)
            log.debug("GET {}", url)

            def request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .GET()
                .build()

            def response = httpClient.sendAsString(request)

            if (response.statusCode() == 404) {
                log.debug("No files found for {}", url)
                return null
            }

            if (response.statusCode() >= 400) {
                log.debug("Error getting files {}", url)
                def message = response.body() ?: ''
                message = "${message} - HTTP ${response.statusCode()}"
                throw new IOException("Failed to get file in '${path}': ${message}")
            }

            def responseJson = jsonSlurper.parseText(response.body()) as Map
            def files = responseJson.objects as List<Map>
            if (files){
                files = files.stream().filter {it.type == 'FOLDER' ? it.name == fileName + '/' : it.name == fileName  } .toList()
            }
            if (files.isEmpty())
                return null
            if (files.size() > 1) {
                throw new IOException("Two items with the same name $fileName")
            }
            Map item = files.get(0)
            return new DataLinkItem(DataLinkItemType.valueOf(item.type as String), item.name as String, item.size as Integer, item.mimeType as String)

        } catch (InterruptedException e) {
            Thread.currentThread().interrupt()
            throw new IOException("List files interrupted", e)
        }
    }

    /**
     * Delete a file from a data-link
     *
     * @param httpClient HTTP client to use
     * @param endpoint Platform API endpoint
     * @param dataLink Data-link containing the file
     * @param filePath Path of file to delete
     * @param workspaceId Workspace ID (optional)
     * @throws IOException if delete fails
     */
    static void deleteFile(HxClient httpClient, String endpoint, DataLink dataLink,
                          String filePath, String workspaceId) throws IOException {
        try {
            final queryParams = [:] as Map<String, String>
            if (workspaceId) {
                queryParams.workspaceId = workspaceId
            }
            if (dataLink.getCredentialsId()) {
                queryParams.credentialsId = dataLink.getCredentialsId()
            }

            def url = buildDataLinkUrl(endpoint, dataLink.getId(), '/content', queryParams)

            def deleteRequest = [
                files: [filePath]
            ]

            def request = HttpRequest.newBuilder()
                .uri(URI.create(url))
                .header('Content-Type', 'application/json')
                .method('DELETE', HttpRequest.BodyPublishers.ofString(new JsonBuilder(deleteRequest).toString()))
                .build()

            def response = httpClient.sendAsString(request)

             if (response.statusCode() >= 400) {
                log.debug("Failed to delete file {}", url)
                def message = response.body() ?: ''
                message = "${message} - HTTP ${response.statusCode()}"
                throw new IOException("Failed to delete file ${filePath}: ${message}")
            }

            def responseJson = jsonSlurper.parseText(response.body()) as Map
            if (responseJson.deletionFailures) {
                def failures = responseJson.deletionFailures as List<Map>
                if (failures) {
                    throw new IOException("Detected deletion failures: ${parseFailures(failures)}")
                }
            }

        } catch (InterruptedException e) {
            Thread.currentThread().interrupt()
            throw new IOException("Delete interrupted", e)
        }
    }

    private static String parseFailures(List<Map> failures) {
        def message = []
        for (Map failure : failures ){
            final error = failure.errorMessage as String
            final item = failure.dataLinkItem as Map
            message.add("${item.type} ${item.name}: ${error}")
        }
        return message.join(',')
    }

    @Canonical
    static class DataLinkItem {
        DataLinkItemType type
        String name
        Integer size
        String mimeType
    }

    enum DataLinkItemType {
        FOLDER,
        FILE,
    }
}
