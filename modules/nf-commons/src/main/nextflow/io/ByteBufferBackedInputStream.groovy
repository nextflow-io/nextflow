/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
package nextflow.io

import java.nio.ByteBuffer

import groovy.transform.CompileStatic

/**
 * An {@code InputStream} adaptor which reads data from a {@code ByteBuffer}
 *
 * See http://stackoverflow.com/a/6603018/395921
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ByteBufferBackedInputStream extends InputStream {

    ByteBuffer buf;

    ByteBufferBackedInputStream(ByteBuffer buf) {
        this.buf = buf;
    }

    int read() throws IOException {
        if (!buf.hasRemaining()) {
            return -1;
        }
        return buf.get() & 0xFF;
    }

    int read(byte[] bytes, int off, int len) throws IOException {
        if (!buf.hasRemaining()) {
            return -1;
        }

        len = Math.min(len, buf.remaining());
        buf.get(bytes, off, len);
        return len;
    }
}