/*
 * Copyright 2013-2026, Seqera Labs
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

package nextflow.cloud.aws.nio.util;

/**
 * @author Paolo Di Tommaso paolo.ditommaso@gmail.com
 */

import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;

/**
 * An {@code InputStream} adaptor which reads data from a {@code ByteBuffer}
 *
 * See http://stackoverflow.com/a/6603018/395921
 *
 * @author Paolo Di Tommaso paolo.ditommaso@gmail.com
 */
public class ByteBufferInputStream extends InputStream {

    ByteBuffer buf;

    /**
     * Position to rewind to on {@link #reset()}: the mark when one was set, otherwise the position
     * the stream started at.
     */
    private int mark;

    public ByteBufferInputStream(ByteBuffer buf) {
        this.buf = buf;
        this.mark = buf.position();
    }

    /**
     * Mark/reset IS supported: the content is a buffer already in memory, so rewinding is a position
     * assignment. It matters because the AWS SDK asks a {@code RequestBody} for a fresh stream on
     * every attempt: over a stream that reports {@code markSupported() == false} it builds a
     * one-shot provider, and a RETRY then fails with "Content input stream does not support
     * mark/reset, and was already read once" instead of re-sending. Both {@code S3OutputStream}
     * upload paths pass this stream -- single-part {@code putObject} and multipart
     * {@code uploadPart} -- so without mark/reset one transient 5xx on any driver-side S3 write
     * aborted the operation instead of being retried.
     */
    @Override
    public boolean markSupported() {
        return true;
    }

    @Override
    public synchronized void mark(int readlimit) {
        this.mark = buf.position();
    }

    @Override
    public synchronized void reset() {
        buf.position(mark);
    }

    @Override
    public int available() {
        return buf.remaining();
    }

    public int read() throws IOException {
        if (!buf.hasRemaining()) {
            return -1;
        }
        return buf.get() & 0xFF;
    }

    public int read(byte[] bytes, int off, int len) throws IOException {
        if (!buf.hasRemaining()) {
            return -1;
        }

        len = Math.min(len, buf.remaining());
        buf.get(bytes, off, len);
        return len;
    }
}
