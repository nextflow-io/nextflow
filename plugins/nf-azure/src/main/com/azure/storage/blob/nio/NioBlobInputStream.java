// Copyright (c) Microsoft Corporation. All rights reserved.
// Licensed under the MIT License.

package com.azure.storage.blob.nio;

import com.azure.core.util.logging.ClientLogger;
import com.azure.storage.blob.specialized.BlobInputStream;

import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Path;

/**
 * Provides an InputStream to read a file stored as an Azure Blob.
 */
public final class NioBlobInputStream extends InputStream {
    private final ClientLogger logger = new ClientLogger(NioBlobInputStream.class);

    private final BlobInputStream blobInputStream;
    private final Path path;

    NioBlobInputStream(BlobInputStream blobInputStream, Path path) {
        this.blobInputStream = blobInputStream;
        this.path = path;
    }

    /**
     * Returns an estimate of the number of bytes that can be read (or skipped over) from this input stream without
     * blocking by the next invocation of a method for this input stream. The next invocation might be the same thread
     * or another thread. A single read or skip of this many bytes will not block, but may read or skip fewer bytes.
     *
     * @return An <code>int</code> which represents an estimate of the number of bytes that can be read (or skipped
     * over) from this input stream without blocking, or 0 when it reaches the end of the input stream.
     */
    @Override
    public synchronized int available() throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        return this.blobInputStream.available();
    }

    /**
     * Closes this input stream and releases any system resources associated with the stream.
     */
    @Override
    public synchronized void close() throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        this.blobInputStream.close();
    }

    /**
     * Marks the current position in this input stream. A subsequent call to the reset method repositions this stream at
     * the last marked position so that subsequent reads re-read the same bytes.
     *
     * @param readlimit An <code>int</code> which represents the maximum limit of bytes that can be read before the mark
     * position becomes invalid.
     */
    @Override
    public synchronized void mark(final int readlimit) {
        this.blobInputStream.mark(readlimit);
    }

    /**
     * Tests if this input stream supports the mark and reset methods.
     *
     * @return Returns {@code true}
     */
    @Override
    public boolean markSupported() {
        return this.blobInputStream.markSupported();
    }

    /**
     * Reads the next byte of data from the input stream. The value byte is returned as an int in the range 0 to 255. If
     * no byte is available because the end of the stream has been reached, the value -1 is returned. This method blocks
     * until input data is available, the end of the stream is detected, or an exception is thrown.
     *
     * @return An <code>int</code> which represents the total number of bytes read into the buffer, or -1 if there is no
     * more data because the end of the stream has been reached.
     * @throws IOException If an I/O error occurs.
     */
    @Override
    public int read() throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        try {
            return this.blobInputStream.read();
            /*
            BlobInputStream only throws RuntimeException, and it doesn't preserve the cause, it only takes the message,
            so we can't do any better than re-wrapping it in an IOException.
             */
        } catch (RuntimeException e) {
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    /**
     * Reads some number of bytes from the input stream and stores them into the buffer array <code>b</code>. The number
     * of bytes actually read is returned as an integer. This method blocks until input data is available, end of file
     * is detected, or an exception is thrown. If the length of <code>b</code> is zero, then no bytes are read and 0 is
     * returned; otherwise, there is an attempt to read at least one byte. If no byte is available because the stream is
     * at the end of the file, the value -1 is returned; otherwise, at least one byte is read and stored into
     * <code>b</code>.
     *
     * The first byte read is stored into element <code>b[0]</code>, the next one into <code>b[1]</code>, and so on. The
     * number of bytes read is, at most, equal to the length of <code>b</code>. Let <code>k</code> be the number of
     * bytes actually read; these bytes will be stored in elements <code>b[0]</code> through <code>b[k-1]</code>,
     * leaving elements <code>b[k]</code> through
     * <code>b[b.length-1]</code> unaffected.
     *
     * The <code>read(b)</code> method for class {@link InputStream} has the same effect as:
     *
     * <code>read(b, 0, b.length)</code>
     *
     * @param b A <code>byte</code> array which represents the buffer into which the data is read.
     * @throws IOException If the first byte cannot be read for any reason other than the end of the file, if the input
     * stream has been closed, or if some other I/O error occurs.
     * @throws NullPointerException If the <code>byte</code> array <code>b</code> is null.
     */
    @Override
    public int read(final byte[] b) throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        try {
            return this.blobInputStream.read(b);
        } catch (RuntimeException e) {
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    /**
     * Reads up to <code>len</code> bytes of data from the input stream into an array of bytes. An attempt is made to
     * read as many as <code>len</code> bytes, but a smaller number may be read. The number of bytes actually read is
     * returned as an integer. This method blocks until input data is available, end of file is detected, or an
     * exception is thrown.
     *
     * If <code>len</code> is zero, then no bytes are read and 0 is returned; otherwise, there is an attempt to read at
     * least one byte. If no byte is available because the stream is at end of file, the value -1 is returned;
     * otherwise, at least one byte is read and stored into <code>b</code>.
     *
     * The first byte read is stored into element <code>b[off]</code>, the next one into <code>b[off+1]</code>, and so
     * on. The number of bytes read is, at most, equal to <code>len</code>. Let <code>k</code> be the number of bytes
     * actually read; these bytes will be stored in elements <code>b[off]</code> through <code>b[off+k-1]</code>,
     * leaving elements <code>b[off+k]</code> through
     * <code>b[off+len-1]</code> unaffected.
     *
     * In every case, elements <code>b[0]</code> through <code>b[off]</code> and elements <code>b[off+len]</code>
     * through <code>b[b.length-1]</code> are unaffected.
     *
     * @param b A <code>byte</code> array which represents the buffer into which the data is read.
     * @param off An <code>int</code> which represents the start offset in the <code>byte</code> array at which the data
     * is written.
     * @param len An <code>int</code> which represents the maximum number of bytes to read.
     * @return An <code>int</code> which represents the total number of bytes read into the buffer, or -1 if there is no
     * more data because the end of the stream has been reached.
     * @throws IOException If the first byte cannot be read for any reason other than end of file, or if the input
     * stream has been closed, or if some other I/O error occurs.
     * @throws NullPointerException If the <code>byte</code> array <code>b</code> is null.
     * @throws IndexOutOfBoundsException If <code>off</code> is negative, <code>len</code> is negative, or
     * <code>len</code> is greater than
     * <code>b.length - off</code>.
     */
    @Override
    public int read(final byte[] b, final int off, final int len) throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        if (off < 0 || len < 0 || len > b.length - off) {
            throw logger.logExceptionAsError(new IndexOutOfBoundsException());
        }
        try {
            return this.blobInputStream.read(b, off, len);
        } catch (RuntimeException e) {
            throw LoggingUtility.logError(logger, new IOException(e));
        }
    }

    /**
     * Repositions this stream to the position at the time the mark method was last called on this input stream. Note
     * repositioning the blob read stream will disable blob MD5 checking.
     *
     * @throws IOException If this stream has not been marked or if the mark has been invalidated.
     */
    @Override
    public synchronized void reset() throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        try {
            this.blobInputStream.reset();
        } catch (RuntimeException e) {
            if (e.getMessage().equals("Stream mark expired.")) {
                throw LoggingUtility.logError(logger, new IOException(e));
            }
            throw LoggingUtility.logError(logger, e);
        }
    }

    /**
     * Skips over and discards n bytes of data from this input stream. The skip method may, for a variety of reasons,
     * end up skipping over some smaller number of bytes, possibly 0. This may result from any of a number of
     * conditions; reaching end of file before n bytes have been skipped is only one possibility. The actual number of
     * bytes skipped is returned. If n is negative, no bytes are skipped.
     *
     * Note repositioning the blob read stream will disable blob MD5 checking.
     *
     * @param n A <code>long</code> which represents the number of bytes to skip.
     */
    @Override
    public synchronized long skip(final long n) throws IOException {
        AzurePath.ensureFileSystemOpen(path);
        return this.blobInputStream.skip(n);
    }
}
