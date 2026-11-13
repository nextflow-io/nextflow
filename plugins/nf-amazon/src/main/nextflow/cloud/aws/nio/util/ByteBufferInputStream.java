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
 * <p>Not thread-safe: the underlying {@code ByteBuffer} position is shared mutable state, so a
 * single instance must be used by one thread at a time. The {@code synchronized} on
 * {@link #mark(int)} / {@link #reset()} follows the {@link InputStream} declaration and does not
 * make the class concurrent.
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
     * Mark/reset is supported: the content is a buffer already in memory, so rewinding is a
     * position assignment. The AWS SDK builds a one-shot content provider over any stream that
     * reports {@code false} here, which makes the request body unusable on a retry.
     *
     * <p>{@code readlimit} is deliberately ignored. The SDK marks with {@code 1 << 17} (128 KiB)
     * while a multipart part is at least 5 MB, so honouring the limit would refuse the reset this
     * override exists to allow.
     *
     * <p>{@link #reset()} never throws, which {@link InputStream} permits: with no prior
     * {@link #mark(int)} it rewinds to the position the stream started at.
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
