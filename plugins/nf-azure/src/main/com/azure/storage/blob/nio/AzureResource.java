// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

import com.azure.core.util.logging.ClientLogger;
import com.azure.storage.blob.BlobContainerClient;
import com.azure.storage.blob.BlobClient;
import com.azure.storage.blob.BlobContainerClientBuilder;
import com.azure.storage.blob.models.BlobHttpHeaders;
import com.azure.storage.blob.models.BlobItem;
import com.azure.storage.blob.models.BlobListDetails;
import com.azure.storage.blob.models.BlobRequestConditions;
import com.azure.storage.blob.models.BlobStorageException;
import com.azure.storage.blob.models.ListBlobsOptions;
import com.azure.storage.common.implementation.Constants;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.attribute.FileAttribute;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Objects;

/**
 * This type is meant to be a logical grouping of operations and data associated with an azure resource. It is NOT
 * intended to serve as a local cache for any data related to remote resources. It is agnostic to whether the resource
 * is a directory or a file and will not perform any validation of the resource type, though root directories are not
 * supported as they are backed by containers and do not support many file system apis.
 *
 * It also serves as the interface to Storage clients. Any operation that needs to use a client should first build an
 * AzureResource using a path and then use the getter to access the client.
 */
final class AzureResource {
    private final ClientLogger logger = new ClientLogger(AzureResource.class);

    static final String DIR_METADATA_MARKER = Constants.HeaderConstants.DIRECTORY_METADATA_KEY;

    private final AzurePath path;
    private final BlobClient blobClient;

    // The following are not kept consistent with the service. They are only held here between parsing and putting.
    private BlobHttpHeaders blobHeaders;
    private Map<String, String> blobMetadata;

    AzureResource(Path path) throws IOException {
        Objects.requireNonNull(path, "path");
        this.path = validatePathInstanceType(path);
        this.validateNotRoot();
        this.blobClient = this.path.toBlobClient();
    }

    /**
     * Checks for the existence of the parent of the given path. We do not check for the actual marker blob as parents
     * need only weakly exist.
     *
     * If the parent is a root (container), it will be assumed to exist, so it must be validated elsewhere that the
     * container is a legitimate root within this file system.
     */
    boolean checkParentDirectoryExists() throws IOException {
        /*
        If the parent is just the root (or null, which means the parent is implicitly the default directory which is a
        root), that means we are checking a container, which is always considered to exist. Otherwise, perform normal
        existence check.
         */
        Path parent = this.path.getParent();
        return (parent == null || parent.equals(path.getRoot()))
            || new AzureResource(this.path.getParent()).checkDirectoryExists();
    }

    /**
     * Checks whether a directory exists by either being empty or having children.
     */
    boolean checkDirectoryExists() throws IOException {
        DirectoryStatus dirStatus = this.checkDirStatus();
        return dirStatus.equals(DirectoryStatus.EMPTY) || dirStatus.equals(DirectoryStatus.NOT_EMPTY);
    }

    /**
     * This method will check if a directory is extant and/or empty and accommodates virtual directories. This method
     * will not check the status of root directories.
     */
    DirectoryStatus checkDirStatus() throws IOException {
        if (this.blobClient == null) {
            throw LoggingUtility.logError(logger, new IllegalArgumentException("The blob client was null."));
        }
        BlobContainerClient containerClient = this.getContainerClient();

        // Two blobs will give us all the info we need (see below).
        ListBlobsOptions listOptions = new ListBlobsOptions().setMaxResultsPerPage(2)
            .setPrefix(this.blobClient.getBlobName())
            .setDetails(new BlobListDetails().setRetrieveMetadata(true));

        /*
        Do a list on prefix.
        Zero elements means no virtual dir. Does not exist.
        One element that matches this dir means empty.
        One element that doesn't match this dir or more than one element. Not empty.
        One element that matches the name but does not have a directory marker means the resource is not a directory.

        Note that blob names that match the prefix exactly are returned in listing operations.
         */
        try {
            Iterator<BlobItem> blobIterator = containerClient.listBlobsByHierarchy(AzureFileSystem.PATH_SEPARATOR,
                listOptions, null).iterator();
            if (!blobIterator.hasNext()) { // Nothing there
                return DirectoryStatus.DOES_NOT_EXIST;
            } else {
                BlobItem item = blobIterator.next();
                if (blobIterator.hasNext()) { // More than one item with dir path as prefix. Must be a dir.
                    return DirectoryStatus.NOT_EMPTY;
                }
                if (!item.getName().equals(this.blobClient.getBlobName())) {
                    /*
                    Names do not match. Must be a virtual dir with one item. e.g. blob with name "foo/bar" means dir
                    "foo" exists.
                     */
                    return DirectoryStatus.NOT_EMPTY;
                }
                if (item.getMetadata() != null && item.getMetadata().containsKey(DIR_METADATA_MARKER)) {
                    return DirectoryStatus.EMPTY; // Metadata marker.
                }
                return DirectoryStatus.NOT_A_DIRECTORY; // There is a file (not a directory) at this location.
            }
        } catch (BlobStorageException e) {
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    /**
     * Creates the actual directory marker. This method should only be used when any necessary checks for proper
     * conditions of directory creation (e.g. parent existence) have already been performed. Otherwise,
     * {@link AzureFileSystemProvider#createDirectory(Path, FileAttribute[])} should be preferred.
     *
     * @param requestConditions Any necessary request conditions to pass when creating the directory blob.
     */
    void putDirectoryBlob(BlobRequestConditions requestConditions) {
        this.blobClient.getBlockBlobClient().commitBlockListWithResponse(Collections.emptyList(), this.blobHeaders,
            this.prepareMetadataForDirectory(), null, requestConditions, null, null);
    }

    /*
    Note that this will remove the properties from the list of attributes as it finds them.
     */
    private void extractHttpHeaders(List<FileAttribute<?>> fileAttributes) {
        BlobHttpHeaders headers = new BlobHttpHeaders();
        for (Iterator<FileAttribute<?>> it = fileAttributes.iterator(); it.hasNext();) {
            FileAttribute<?> attr = it.next();
            boolean propertyFound = true;
            switch (attr.name()) {
                case AzureFileSystemProvider.CONTENT_TYPE:
                    headers.setContentType(attr.value().toString());
                    break;
                case AzureFileSystemProvider.CONTENT_LANGUAGE:
                    headers.setContentLanguage(attr.value().toString());
                    break;
                case AzureFileSystemProvider.CONTENT_DISPOSITION:
                    headers.setContentDisposition(attr.value().toString());
                    break;
                case AzureFileSystemProvider.CONTENT_ENCODING:
                    headers.setContentEncoding(attr.value().toString());
                    break;
                case AzureFileSystemProvider.CONTENT_MD5:
                    if ((attr.value() instanceof byte[])) {
                        headers.setContentMd5((byte[]) attr.value());
                    } else {
                        throw LoggingUtility.logError(logger,
                            new UnsupportedOperationException("Content-MD5 attribute must be a byte[]"));
                    }
                    break;
                case AzureFileSystemProvider.CACHE_CONTROL:
                    headers.setCacheControl(attr.value().toString());
                    break;
                default:
                    propertyFound = false;
                    break;
            }

            if (propertyFound) {
                it.remove();
            }
        }

        this.blobHeaders = headers;
    }

    /**
     * Note this should only be used after the headers have been extracted.
     *
     * @param fileAttributes The attributes to convert to metadata.
     */
    private void convertAttributesToMetadata(List<FileAttribute<?>> fileAttributes) {
        Map<String, String> metadata = new HashMap<>();
        for (FileAttribute<?> attr : fileAttributes) {
            metadata.put(attr.name(), attr.value().toString());
        }

        // If no attributes are set, return null so existing metadata is not cleared.
        this.blobMetadata = metadata.isEmpty() ? null : metadata;
    }

    private void validateNotRoot() {
        if (this.path.isRoot()) {
            throw LoggingUtility.logError(logger, new IllegalArgumentException(
                "Root directory not supported. Path: " + this.path.toString()));
        }
    }

    private AzurePath validatePathInstanceType(Path path) {
        if (!(path instanceof AzurePath)) {
            throw LoggingUtility.logError(logger, new IllegalArgumentException("This provider cannot operate on "
                + "subtypes of Path other than AzurePath"));
        }
        return (AzurePath) path;
    }

    BlobContainerClient getContainerClient() {
        return new BlobContainerClientBuilder().endpoint(this.blobClient.getBlobUrl())
            .pipeline(this.blobClient.getHttpPipeline())
            .buildClient();
    }

    AzureResource setFileAttributes(List<FileAttribute<?>> attributes) {
        attributes = new ArrayList<>(attributes); // To ensure removing header values from the list is supported.
        extractHttpHeaders(attributes);
        convertAttributesToMetadata(attributes);

        return this;
    }

    AzurePath getPath() {
        return this.path;
    }

    BlobClient getBlobClient() {
        return this.blobClient;
    }

    private Map<String, String> prepareMetadataForDirectory() {
        if (this.blobMetadata == null) {
            this.blobMetadata = new HashMap<>();
        }
        this.blobMetadata.put(DIR_METADATA_MARKER, "true");
        return this.blobMetadata;
    }
}
