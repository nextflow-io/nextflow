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

import java.time.Duration
import java.time.temporal.ChronoUnit

import dev.failsafe.ExecutionContext
import dev.failsafe.Failsafe
import dev.failsafe.RetryPolicy
import dev.failsafe.function.ContextualSupplier
import groovy.util.logging.Slf4j
import spock.lang.Ignore
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
class ChunkedInputStreamTest extends Specification {

    def 'should read the chunks' () {
        given:
        def chunk0 = new String('Hello world\n').bytes
        def chunk1 = new String('Hola mundo\n').bytes
        def chunk2 = new String('Ciao mondo\n').bytes
        and:
        def len = chunk0.length+chunk1.length+chunk2.length
        def stream = new ChunkedInputStream(len)

        when:
        stream.add(ChunkBuffer.wrap(chunk0).withIndex(0))
        stream.add(ChunkBuffer.wrap(chunk1).withIndex(1))
        stream.add(ChunkBuffer.wrap(chunk2).withIndex(2))

        then:
        stream.text == '''\
            Hello world
            Hola mundo
            Ciao mondo
            '''.stripIndent()
    }

    def 'should read the chunks async' () {
        given:
        def chunk0 = new String('Hello world\n').bytes
        def chunk1 = new String('Hola mundo\n').bytes
        def chunk2 = new String('Ciao mondo\n').bytes
        and:
        def len = chunk0.length+chunk1.length+chunk2.length
        def stream = new ChunkedInputStream(len)

        when:
        Thread.start { sleep 100; stream.add(ChunkBuffer.wrap(chunk0).withIndex(0)) }
        Thread.start { sleep 200; stream.add(ChunkBuffer.wrap(chunk1).withIndex(1)) }
        Thread.start { sleep 300; stream.add(ChunkBuffer.wrap(chunk2).withIndex(2)) }

        then:
        stream.text == '''\
            Hello world
            Hola mundo
            Ciao mondo
            '''.stripIndent()
    }

    def 'should read empty string' () {
        given:
        def stream = new ChunkedInputStream(0)
        expect:
        stream.text == ''
    }

    def 'should read throw an exception' () {
        given:
        def chunk0 = new String('Hello world\n').bytes
        def chunk1 = new String('Hola mundo\n').bytes
        def chunk2 = new String('Ciao mondo\n').bytes
        and:
        def len = chunk0.length+chunk1.length+chunk2.length
        def stream = new ChunkedInputStream(len)

        when:
        Thread.start { sleep 100; stream.add(ChunkBuffer.wrap(chunk0)) }
        Thread.start { sleep 200; stream.throwError(new IOException("Something break")) }
        and:
        println stream.text

        then:
        thrown(IOException)
    }

    def 'should read the stream ad give bakc the chunks' () {
        given:
        def STR = "hello world!"
        def BYTES = STR.bytes
        def CHUNK_SIZE = BYTES.length  +2
        def POOL_CAPACITY = 10
        def TIMES = 10
        def buffers = new ChunkBufferFactory(CHUNK_SIZE, POOL_CAPACITY)
        and:
        def LEN = BYTES.length * TIMES;
        and:
        def executor = PriorityThreadPool.create('foo', 10, 1000)

        when:
        def stream = new ChunkedInputStream(LEN)
        and:
        TIMES.times { index ->
            executor.submit( new PriorityThreadPool.PriorityRunnable(index) {
                @Override
                void run() {
                    def chunk = buffers.create(index)
                    chunk.fill( new ByteArrayInputStream(BYTES) )
                    chunk.makeReadable()
                    stream.add(chunk)
                }
            })
        }

        then:
        stream.text == STR * TIMES
        and:
        buffers.getPoolSize() == POOL_CAPACITY

        cleanup:
        executor.shutdownNow()
    }

    @Ignore
    def 'test failsafe' () {
        given:
        RetryPolicy<Object> retryPolicy = RetryPolicy.builder()
                .handle(RuntimeException.class)
                .withDelay(Duration.ofSeconds(1))
                .withMaxDuration(Duration.of(60, ChronoUnit.SECONDS))
                .withBackoff(1, 30, ChronoUnit.SECONDS)
                .withMaxRetries(10)
                .onFailedAttempt(e -> log.error("Connection attempt failed - cause: ${e.getLastFailure()}"))
                .onRetry(e -> log.warn("Failure #{}. Retrying.", e.getAttemptCount()))
                .build();

        when:
        def work = { ExecutionContext it ->
            log.debug "try num ${it.getAttemptCount()}"
            throw new RuntimeException("Break ${it.getAttemptCount()}")
            } as ContextualSupplier
        def result = Failsafe.with(retryPolicy).get( work )
        then:
        result == 'Hello'
    }

}
