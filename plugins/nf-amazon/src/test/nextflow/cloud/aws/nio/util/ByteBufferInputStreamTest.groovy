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

package nextflow.cloud.aws.nio.util

import java.nio.ByteBuffer

import software.amazon.awssdk.core.sync.RequestBody
import spock.lang.Specification

class ByteBufferInputStreamTest extends Specification {

    def 'should read the buffer content'() {
        given:
        def stream = new ByteBufferInputStream(ByteBuffer.wrap('hello'.bytes))
        expect:
        stream.bytes == 'hello'.bytes
    }

    /**
     * The SDK asks a {@code RequestBody} for a fresh stream on every attempt, so a body that can be
     * read only once cannot be RETRIED: the second attempt fails with IllegalStateException instead
     * of re-sending.
     */
    def 'a request body built from the buffer can be re-read, so the SDK can retry'() {
        given:
        def payload = 'index-bytes-for-the-cache'.bytes
        def body = RequestBody.fromInputStream(new ByteBufferInputStream(ByteBuffer.wrap(payload)), payload.length)

        when: 'the SDK asks for the stream twice, as it does when the first attempt failed'
        def first = body.contentStreamProvider().newStream().bytes
        def second = body.contentStreamProvider().newStream().bytes

        then: 'both attempts see the whole payload'
        first == payload
        second == payload
    }

    def 'a partially consumed stream rewinds to the mark'() {
        given: 'the shape of a failed attempt -- some bytes went out before the error'
        def payload = '0123456789'.bytes
        def stream = new ByteBufferInputStream(ByteBuffer.wrap(payload))

        expect:
        stream.markSupported()

        when:
        stream.mark(payload.length)
        def partial = new byte[4]
        stream.read(partial, 0, 4)
        stream.reset()

        then: 'the retry re-sends the whole body, not the tail'
        partial == '0123'.bytes
        stream.bytes == payload
    }

    def 'reset without an explicit mark rewinds to the start'() {
        given: 'the SDK may reset a provider it never marked'
        def stream = new ByteBufferInputStream(ByteBuffer.wrap('abc'.bytes))

        when:
        stream.read()
        stream.reset()

        then:
        stream.bytes == 'abc'.bytes
    }

    def 'available reports the unread remainder'() {
        given:
        def stream = new ByteBufferInputStream(ByteBuffer.wrap('abcde'.bytes))
        expect:
        stream.available() == 5
        when:
        stream.read()
        then:
        stream.available() == 4
    }
}
