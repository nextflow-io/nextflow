// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

import com.azure.core.util.logging.ClientLogger;
import com.azure.storage.blob.models.BlobProperties;
import com.azure.storage.blob.models.BlobStorageException;

import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.attribute.BasicFileAttributes;
import java.nio.file.attribute.FileAttribute;
import java.nio.file.attribute.FileTime;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 * Provides support for basic file attributes.
 * <p>
 * The properties available on this type are a strict subset of {@link AzureBlobFileAttributes}, and the two types have
 * the same network behavior. Therefore, while this type is offered for compliance with the NIO spec,
 * {@link AzureBlobFileAttributes} is generally preferred.
 * <p>
 * Some attributes are not supported. Refer to the javadocs on each method for more information.
 */
public final class AzureBasicFileAttributes implements BasicFileAttributes {
    private final ClientLogger logger = new ClientLogger(AzureBasicFileAttributes.class);

    // For verifying parameters on FileSystemProvider.readAttributes
    static final Set<String> ATTRIBUTE_STRINGS;
    static {
        Set<String> set = new HashSet<>();
        set.add("lastModifiedTime");
        set.add("isRegularFile");
        set.add("isDirectory");
        set.add("isSymbolicLink");
        set.add("isOther");
        set.add("size");
        set.add("creationTime");
        ATTRIBUTE_STRINGS = Collections.unmodifiableSet(set);
    }

    private final BlobProperties properties;
    private final AzureResource resource;

    /*
    There are some work-arounds we could do to try to accommodate virtual directories such as making a checkDirStatus
    call before or after getProperties to throw an appropriate error or adding an isVirtualDirectory method. However,
    the former wastes network time only to throw a slightly more specific error when we will throw on 404 anyway. The
    latter introduces virtual directories into the actual code path/api surface. While we are clear in our docs about
    the possible pitfalls of virtual directories, and customers should be aware of it, they shouldn't have to code
    against it. Therefore, we fall back to documenting that reading attributes on a virtual directory will throw.
     */
    AzureBasicFileAttributes(Path path) throws IOException {
        try {
            this.resource = new AzureResource(path);
            this.properties = resource.getBlobClient().getProperties();
        } catch (BlobStorageException e) {
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    /**
     * Returns the time of last modification.
     *
     * @return the time of last modification.
     */
    @Override
    public FileTime lastModifiedTime() {
        return FileTime.from(properties.getLastModified().toInstant());
    }

    /**
     * Returns the time of last modification.
     * <p>
     * Last access time is not supported by the blob service. In this case, it is typical for implementations to return
     * the {@link #lastModifiedTime()}.
     *
     * @return the time of last modification.
     */
    @Override
    public FileTime lastAccessTime() {
        return this.lastModifiedTime();
    }

    /**
     * Returns the creation time. The creation time is the time that the file was created.
     *
     * @return The creation time.
     */
    @Override
    public FileTime creationTime() {
        return FileTime.from(properties.getCreationTime().toInstant());
    }

    /**
     * Tells whether the file is a regular file with opaque content.
     *
     * @return whether the file is a regular file.
     */
    @Override
    public boolean isRegularFile() {
        return !this.properties.getMetadata().getOrDefault(AzureResource.DIR_METADATA_MARKER, "false").equals("true");
    }

    /**
     * Tells whether the file is a directory.
     * <p>
     * Will only return true if the directory is a concrete directory. See
     * {@link AzureFileSystemProvider#createDirectory(Path, FileAttribute[])} for more information on virtual and
     * concrete directories.
     *
     * @return whether the file is a directory
     */
    @Override
    public boolean isDirectory() {
        return !this.isRegularFile();
    }

    /**
     * Tells whether the file is a symbolic link.
     *
     * @return false. Symbolic links are not supported.
     */
    @Override
    public boolean isSymbolicLink() {
        return false;
    }

    /**
     * Tells whether the file is something other than a regular file, directory, or symbolic link.
     *
     * @return false. No other object types are supported.
     */
    @Override
    public boolean isOther() {
        return false;
    }

    /**
     * Returns the size of the file (in bytes).
     *
     * @return the size of the file
     */
    @Override
    public long size() {
        return properties.getBlobSize();
    }

    /**
     * Returns the url of the resource.
     *
     * @return The file key, which is the url.
     */
    @Override
    public Object fileKey() {
        return this.resource.getBlobClient().getBlobUrl();
    }
}
