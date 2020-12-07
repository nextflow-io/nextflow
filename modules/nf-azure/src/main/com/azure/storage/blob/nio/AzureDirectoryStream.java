// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

import com.azure.core.util.logging.ClientLogger;
import com.azure.storage.blob.BlobContainerClient;
import com.azure.storage.blob.models.BlobItem;
import com.azure.storage.blob.models.BlobListDetails;
import com.azure.storage.blob.models.ListBlobsOptions;

import java.io.IOException;
import java.nio.file.DirectoryIteratorException;
import java.nio.file.DirectoryStream;
import java.nio.file.Path;
import java.util.HashSet;
import java.util.Iterator;
import java.util.NoSuchElementException;
import java.util.Set;

/**
 * A type for iterating over the contents of a directory.
 *
 * This type is asynchronously closeable, i.e. closing the stream from any thread will cause the stream to stop
 * returning elements at that point.
 *
 * {@inheritDoc}
 */
public final class AzureDirectoryStream implements DirectoryStream<Path> {
    private final ClientLogger logger = new ClientLogger(AzureDirectoryStream.class);

    private final AzurePath path;
    private final DirectoryStream.Filter<? super Path> filter;
    private boolean iteratorRequested = false;
    private final AzureDirectoryIterator iterator;
    boolean closed = false;

    AzureDirectoryStream(AzurePath path, DirectoryStream.Filter<? super Path> filter) throws IOException {
        this.path = path;
        this.filter = filter;
        this.iterator = new AzureDirectoryIterator(this, this.path, this.filter);
    }

    @Override
    public Iterator<Path> iterator() {
        if (this.iteratorRequested) {
            throw LoggingUtility.logError(logger,
                new IllegalStateException("Only one iterator may be requested from a given directory stream"));
        }
        this.iteratorRequested = true;
        return this.iterator;
    }

    @Override
    public void close() throws IOException {
        this.closed = true;
    }

    private static class AzureDirectoryIterator implements Iterator<Path> {
        private final ClientLogger logger = new ClientLogger(AzureDirectoryIterator.class);

        private final AzureDirectoryStream parentStream;
        private final DirectoryStream.Filter<? super Path> filter;
        private final Iterator<BlobItem> blobIterator;
        private final AzurePath path;
        private final Path withoutRoot;
        private Path bufferedNext = null;
        private final Set<String> directoryPaths;

        AzureDirectoryIterator(AzureDirectoryStream parentStream, AzurePath path,
            DirectoryStream.Filter<? super Path> filter) throws IOException {
            this.parentStream = parentStream;
            this.filter = filter;
            this.path = path;

            /*
            Resolving two paths requires that either both have a root or neither does. Because the paths returned from
            listing will never have a root, we prepare a copy of the list path without a root for quick resolving later.
             */
            Path root = this.path.getRoot();
            this.withoutRoot = root == null ? this.path : root.relativize(this.path);

            directoryPaths = new HashSet<>();

            BlobContainerClient containerClient;
            ListBlobsOptions listOptions = new ListBlobsOptions()
                .setDetails(new BlobListDetails().setRetrieveMetadata(true));
            if (path.isRoot()) {
                String containerName = path.toString().substring(0, path.toString().length() - 1);
                AzureFileSystem afs = ((AzureFileSystem) path.getFileSystem());
                containerClient = ((AzureFileStore) afs.getFileStore(containerName)).getContainerClient();
            } else {
                AzureResource azureResource = new AzureResource(path);
                listOptions.setPrefix(azureResource.getBlobClient().getBlobName() + AzureFileSystem.PATH_SEPARATOR);
                containerClient = azureResource.getContainerClient();
            }
            this.blobIterator = containerClient
                .listBlobsByHierarchy(AzureFileSystem.PATH_SEPARATOR, listOptions, null).iterator();
        }

        @Override
        public boolean hasNext() {
            try {
                AzurePath.ensureFileSystemOpen(path);
            } catch (IOException e) {
                throw LoggingUtility.logError(logger, new DirectoryIteratorException(e));
            }
            // Closing the parent stream halts iteration.
            if (parentStream.closed) {
                return false;
            }

            // In case a customer calls hasNext multiple times in a row. If we've buffered an element, we have a next.
            if (this.bufferedNext != null) {
                return true;
            }

            /*
            Search for a new element that passes the filter and buffer it when found. If no such element is found,
            return false.
             */
            while (this.blobIterator.hasNext()) {
                BlobItem nextBlob = this.blobIterator.next();
                Path nextPath = getNextListResult(nextBlob);
                try {
                    if (this.filter.accept(nextPath) && isNotDuplicate(nextPath, nextBlob)) {
                        this.bufferedNext = nextPath;
                        return true;
                    }
                } catch (IOException e) {
                    throw LoggingUtility.logError(logger, new DirectoryIteratorException(e));
                }
            }
            return false;
        }

        @Override
        public Path next() {
            if (this.bufferedNext == null) {
                if (!this.hasNext()) { // This will populate bufferedNext in the process.
                    throw LoggingUtility.logError(logger, new NoSuchElementException());
                }
            }
            Path next = this.bufferedNext; // bufferedNext will have been populated by hasNext()
            this.bufferedNext = null;
            return next;
        }

        @Override
        public void remove() {
            throw new UnsupportedOperationException();
        }

        private Path getNextListResult(BlobItem blobItem) {
            /*
            Listing results return the full blob path, and we don't want to duplicate the path we listed off of, so
            we relativize to remove it.
             */
            String blobName = blobItem.getName();
            Path relativeResult = this.withoutRoot.relativize(
                this.path.getFileSystem().getPath(blobName));

            // Resolve the cleaned list result against the original path for the final result.
            return this.path.resolve(relativeResult);
        }

        /*
        If there is a concrete directory with children, a given path will be returned twice: once as the marker blob
        and once as the prefix for its children. We don't want to return the item twice, and we have no guarantees on
        result ordering, so we have to maintain a cache of directory paths we've seen in order to de-dup.
         */
        private boolean isNotDuplicate(Path path, BlobItem blob) {
            /*
            If the blob is not a prefix and the blob does not contain the directory metadata marker, it is a normal blob
            and therefore will not be duplicated.
             */
            if (!(blob.isPrefix() != null && blob.isPrefix())
                && !(blob.getMetadata() != null && blob.getMetadata().containsKey(AzureResource.DIR_METADATA_MARKER))) {
                return true;
            }

            // If the set contains this path, it means we've seen it before and we shouldn't return it again.
            if (this.directoryPaths.contains(path.toString())) {
                return false;
            }

            // We haven't seen this before. Track it and indicate it should be returned.
            this.directoryPaths.add(path.toString());
            return true;
        }
    }
}
