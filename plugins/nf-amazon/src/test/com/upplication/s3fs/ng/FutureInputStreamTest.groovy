/*
 * Copyright 2020-2021, Seqera Labs
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

import java.util.concurrent.Callable
import java.util.concurrent.Executors

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class FutureInputStreamTest extends Specification {

    def 'should read the stream ad give bakc the chunks' () {
        given:
        def STR = "hello world!"
        def BYTES = STR.bytes
        def CHUNK_SIZE = BYTES.length  +2
        def TIMES = 10
        def POOL_CAPACITY = 10
        def buffers = new ChunkBufferFactory(CHUNK_SIZE, POOL_CAPACITY)
        and:
        def executor = Executors.newFixedThreadPool(10)

        when:
        def futures = []
        for( int it : 1..TIMES ) {
            futures << executor.submit( new Callable<ChunkBuffer>() {
                @Override
                ChunkBuffer call() throws Exception {
                    def chunk = buffers.create(0)
                    chunk.fill( new ByteArrayInputStream(BYTES) )
                    chunk.makeReadable()
                    return chunk
                }
            })
        }
        and:
        def stream = new FutureInputStream(futures)

        then:
        stream.text == STR * TIMES
        and:
        buffers.getPoolSize() == POOL_CAPACITY

        cleanup:
        executor.shutdownNow()
    }


}
