// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

import com.azure.core.util.logging.ClientLogger;
import com.azure.storage.blob.BlobContainerClient;

import java.io.IOException;
import java.nio.file.FileStore;
import java.nio.file.attribute.FileAttributeView;
import java.nio.file.attribute.FileStoreAttributeView;
import java.util.Objects;

/**
 * An {@code AzureFileStore} is a {@link FileStore} backed by an Azure Blob Storage container.
 */
public final class AzureFileStore extends FileStore {
    private final ClientLogger logger = new ClientLogger(AzureFileStore.class);

    private static final String AZURE_FILE_STORE_TYPE = "AzureBlobContainer";

    private final AzureFileSystem parentFileSystem;
    private final BlobContainerClient containerClient;


    AzureFileStore(AzureFileSystem parentFileSystem, String containerName) throws IOException {
        // A FileStore should only ever be created by a FileSystem.
        if (Objects.isNull(parentFileSystem)) {
            throw LoggingUtility.logError(logger, new IllegalStateException("AzureFileStore cannot be instantiated "
                + "without a parent FileSystem"));
        }
        this.parentFileSystem = parentFileSystem;
        this.containerClient = this.parentFileSystem.getBlobServiceClient().getBlobContainerClient(containerName);

        try {
            // This also serves as our connection check.
            if (!this.containerClient.exists()) {
                this.containerClient.create();
            }
        } catch (Exception e) {
            throw LoggingUtility.logError(logger, new IOException("There was an error in establishing the existence of "
                + "container: " + containerName, e));
        }
    }

    /**
     * Returns the name of the container that underlies this file store.
     *
     * @return the name of the container that underlies this file store.
     */
    @Override
    public String name() {
        return this.containerClient.getBlobContainerName();
    }

    /**
     * Returns the {@code String "AzureBlobContainer"} to indicate that the file store is backed by a remote blob
     * container in Azure Storage.
     *
     * @return {@code "AzureBlobContainer"}
     */
    @Override
    public String type() {
        return AZURE_FILE_STORE_TYPE;
    }

    /**
     * Always returns false.
     * <p>
     * It may be the case that the authentication method provided to this file system only
     * supports read operations and hence the file store is implicitly read only in this view, but that does not
     * imply the underlying container/file store is inherently read only. Creating/specifying read only file stores
     * is not currently supported.
     *
     * @return false.
     */
    @Override
    public boolean isReadOnly() {
        return false;
    }

    /**
     * Returns the size, in bytes, of the file store.
     * <p>
     * Containers do not limit the amount of data stored. This method will always return max long.
     *
     * @return the size of the file store.
     * @throws IOException If an I/O error occurs.
     */
    @Override
    public long getTotalSpace() throws IOException {
        return Long.MAX_VALUE;
    }

    /**
     * Returns the number of bytes available to this Java virtual machine on the file store.
     * <p>
     * Containers do not limit the amount of data stored. This method will always return max long.
     *
     * @return the number of bytes available on the file store.
     * @throws IOException If an I/O error occurs.
     */
    @Override
    public long getUsableSpace() throws IOException {
        return Long.MAX_VALUE;
    }

    /**
     * Returns the number of unallocated bytes in the file store.
     * <p>
     * Containers do not limit the amount of data stored. This method will always return max long.
     *
     * @return the number of unallocated bytes in the file store.
     * @throws IOException If an I/O error occurs.
     */
    @Override
    public long getUnallocatedSpace() throws IOException {
        return Long.MAX_VALUE;
    }

    /**
     * Tells whether or not this file store supports the file attributes identified by the given file attribute view.
     * <p>
     * All file stores in this file system support the following views:
     * <ul>
     *     <li>{@link java.nio.file.attribute.BasicFileAttributeView}</li>
     *     <li>{@link AzureBasicFileAttributeView}</li>
     *     <li>{@link AzureBlobFileAttributeView}</li>
     * </ul>
     *
     * @param type the file attribute view type
     * @return Whether the file attribute view is supported.
     */
    @Override
    public boolean supportsFileAttributeView(Class<? extends FileAttributeView> type) {
        return AzureFileSystem.SUPPORTED_ATTRIBUTE_VIEWS.containsKey(type);
    }

    /**
     * Tells whether or not this file store supports the file attributes identified by the given file attribute view.
     * <p>
     * All file stores in this file system support the following views:
     * <ul>
     *     <li>{@link java.nio.file.attribute.BasicFileAttributeView}</li>
     *     <li>{@link AzureBasicFileAttributeView}</li>
     *     <li>{@link AzureBlobFileAttributeView}</li>
     * </ul>
     *
     * @param name the name of the file attribute view
     * @return whether the file attribute view is supported.
     */
    @Override
    public boolean supportsFileAttributeView(String name) {
        return AzureFileSystem.SUPPORTED_ATTRIBUTE_VIEWS.containsValue(name);
    }

    /**
     * Returns a FileStoreAttributeView of the given type.
     * <p>
     * This method always returns null as no {@link FileStoreAttributeView} is currently supported.
     *
     * @param aClass a class
     * @return null
     */
    @Override
    public <V extends FileStoreAttributeView> V getFileStoreAttributeView(Class<V> aClass) {
        return null;
    }

    /**
     * Unsupported.
     * <p>
     * This method always throws an {@code UnsupportedOperationException} as no {@link FileStoreAttributeView} is
     * currently supported.
     *
     * @param s a string
     * @return The attribute value.
     * @throws UnsupportedOperationException unsupported
     * @throws IOException never
     */
    @Override
    public Object getAttribute(String s) throws IOException {
        throw LoggingUtility.logError(logger, new UnsupportedOperationException("FileStoreAttributeViews aren't"
            + " supported."));
    }

    BlobContainerClient getContainerClient() {
        return this.containerClient;
    }
}
