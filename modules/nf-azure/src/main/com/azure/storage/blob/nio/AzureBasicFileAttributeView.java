// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.attribute.BasicFileAttributeView;
import java.nio.file.attribute.FileTime;

/**
 * Provides support for basic file attributes.
 * <p>
 * The operations supported by this view and the attributes it reads are a strict subset of
 * {@link AzureBlobFileAttributeView} and has the same network behavior. Therefore, while this type is offered for
 * compliance with the NIO spec, {@link AzureBlobFileAttributeView} is generally preferred.
 * <p>
 * {@link #setTimes(FileTime, FileTime, FileTime)} is not supported.
 */
public final class AzureBasicFileAttributeView implements BasicFileAttributeView {

    static final String NAME = "azureBasic";

    private final Path path;

    AzureBasicFileAttributeView(Path path) {
        this.path = path;
    }

    /**
     * Returns the name of the attribute view: {@code "azureBasic"}
     *
     * @return the name of the attribute view: {@code "azureBasic"}
     */
    @Override
    public String name() {
        return NAME;
    }

    /**
     * Reads the basic file attributes as a bulk operation.
     * <p>
     * All file attributes are read as an atomic operation with respect to other file system operations.
     *
     * @return {@link AzureBasicFileAttributes}
     */
    @Override
    public AzureBasicFileAttributes readAttributes() throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        return new AzureBasicFileAttributes(path);
    }

    /**
     * Unsupported.
     *
     * @param lastModifiedTime the new last modified time, or null to not change the value
     * @param lastAccessTime the last access time, or null to not change the value
     * @param createTime the file's create time, or null to not change the value
     * @throws UnsupportedOperationException Operation not supported.
     * @throws IOException never
     */
    @Override
    public void setTimes(FileTime lastModifiedTime, FileTime lastAccessTime, FileTime createTime) throws IOException {
        throw new UnsupportedOperationException();
    }
}
