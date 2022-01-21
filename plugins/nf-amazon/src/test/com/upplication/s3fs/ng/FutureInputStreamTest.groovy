/*
 * Copyright 2020-2022, Seqera Labs
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

package com.upplication.s3fs.ng


import java.util.concurrent.Executors
import java.util.function.Function

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FutureInputStreamTest extends Specification {

    def 'should read the stream ad give back the chunks' () {
        given:
        def STR = "hello world!"
        def BYTES = STR.bytes
        def CHUNK_SIZE = BYTES.length  +2
        def TIMES = 10
        def CAPACITY = 1
        def buffers = new ChunkBufferFactory(CHUNK_SIZE, CAPACITY)
        and:
        def executor = Executors.newFixedThreadPool(10)

        and:
        def parts = []; TIMES.times { parts.add(it) }
        def Function<Integer,ChunkBuffer> task = {
            def chunk = buffers.create()
            chunk.fill( new ByteArrayInputStream(BYTES) )
            chunk.makeReadable()
            return chunk
        }

        when:
        def itr = new FutureIterator(parts, task, executor, CAPACITY)
        def stream = new FutureInputStream(itr)

        then:
        stream.text == STR * TIMES
        and:
        buffers.getPoolSize() == CAPACITY

        cleanup:
        executor.shutdownNow()
    }


}
