/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.file.http

import groovy.transform.CompileStatic

/**
 * Implements a {@link FilterInputStream} that checks the expected length of bytes have been
 * read when closing the stream or throws an error otherwise
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class FixedInputStream extends FilterInputStream {

    private final long length
    private long bytesRead

    FixedInputStream(InputStream inputStream, long len) {
        super(inputStream)
        this.length = len
    }

    @Override
    int read() throws IOException {
        final result = super.read()
        if( result!=-1 )
            bytesRead++
        return result
    }

    @Override
    int read(byte[] b, int off, int len) throws IOException {
        final result = super.read(b, off, len)
        if( result!=-1 )
            bytesRead += result
        return result
    }

    @Override
    long skip(long n) throws IOException {
        long skipped = super.skip(n)
        bytesRead += skipped
        return skipped
    }

    @Override
    int available() throws IOException {
        super.available()
    }

    @Override
    void close() throws IOException {
        if( bytesRead != length )
            throw new IOException("Read data length does not match expected size - bytes read: ${bytesRead}; expected: ${length}")
        super.close()
    }
}
