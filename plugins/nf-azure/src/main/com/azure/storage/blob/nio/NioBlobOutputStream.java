// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

import com.azure.core.util.logging.ClientLogger;
import com.azure.storage.blob.specialized.BlobOutputStream;

import java.io.IOException;
import java.io.OutputStream;
import java.nio.file.Path;

/**
 * Provides an OutputStream to write to a file stored as an Azure Blob.
 */
public final class NioBlobOutputStream extends OutputStream {
    private final ClientLogger logger = new ClientLogger(NioBlobOutputStream.class);

    private final BlobOutputStream blobOutputStream;
    private final Path path;

    NioBlobOutputStream(BlobOutputStream blobOutputStream, Path path) {
        this.blobOutputStream = blobOutputStream;
        this.path = path;
    }

    @Override
    public synchronized void write(int i) throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        try {
            this.blobOutputStream.write(i);
            /*
            BlobOutputStream only throws RuntimeException, and it doesn't preserve the cause, it only takes the message,
            so we can't do any better than re-wrapping it in an IOException.
             */
        } catch (RuntimeException e) {
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    @Override
    public synchronized void write(byte[] b) throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        try {
            this.blobOutputStream.write(b);
            /*
            BlobOutputStream only throws RuntimeException, and it doesn't preserve the cause, it only takes the message,
            so we can't do any better than re-wrapping it in an IOException.
             */
        } catch (RuntimeException e) {
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    @Override
    public synchronized void write(byte[] b, int off, int len) throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        try {
            this.blobOutputStream.write(b, off, len);
            /*
            BlobOutputStream only throws RuntimeException, and it doesn't preserve the cause, it only takes the message,
            so we can't do any better than re-wrapping it in an IOException.
             */
        } catch (RuntimeException e) {
            if (e instanceof IndexOutOfBoundsException) {
                throw LoggingUtility.logError(logger, e);
            }
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    @Override
    public synchronized void flush() throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        try {
            this.blobOutputStream.flush();
            /*
            BlobOutputStream only throws RuntimeException, and it doesn't preserve the cause, it only takes the message,
            so we can't do any better than re-wrapping it in an IOException.
             */
        } catch (RuntimeException e) {
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    @Override
    public synchronized void close() throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        try {
            this.blobOutputStream.close();
            /*
            BlobOutputStream only throws RuntimeException, and it doesn't preserve the cause, it only takes the message,
            so we can't do any better than re-wrapping it in an IOException.
             */
        } catch (RuntimeException e) {
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }
}
