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

package io.seqera.tower.plugin.fs

import java.nio.ByteBuffer

import spock.lang.Specification

class DatasetInputStreamTest extends Specification {

    def 'should read bytes into buffer'() {
        given:
        def data = 'hello world'.bytes
        def channel = new DatasetInputStream(new ByteArrayInputStream(data))
        def buf = ByteBuffer.allocate(data.length)

        when:
        def n = channel.read(buf)

        then:
        n == data.length
        buf.array() == data
    }

    def 'should advance position after read'() {
        given:
        def data = 'abcdef'.bytes
        def channel = new DatasetInputStream(new ByteArrayInputStream(data))

        when:
        channel.read(ByteBuffer.allocate(3))

        then:
        channel.position() == 3

        when:
        channel.read(ByteBuffer.allocate(3))

        then:
        channel.position() == 6
    }

    def 'should return -1 at end of stream'() {
        given:
        def channel = new DatasetInputStream(new ByteArrayInputStream(new byte[0]))

        when:
        def n = channel.read(ByteBuffer.allocate(8))

        then:
        n == -1
        channel.position() == 0
    }

    def 'should read partial buffer when stream has fewer bytes'() {
        given:
        def data = 'hi'.bytes
        def channel = new DatasetInputStream(new ByteArrayInputStream(data))
        def buf = ByteBuffer.allocate(100)

        when:
        def n = channel.read(buf)

        then:
        n == 2
        channel.position() == 2
    }

    def 'should be open initially and closed after close()'() {
        given:
        def channel = new DatasetInputStream(new ByteArrayInputStream(new byte[0]))

        expect:
        channel.isOpen()

        when:
        channel.close()

        then:
        !channel.isOpen()
    }

    def 'should close underlying stream on close()'() {
        given:
        def stream = Mock(InputStream)
        def channel = new DatasetInputStream(stream)

        when:
        channel.close()

        then:
        1 * stream.close()
        !channel.isOpen()
    }

    def 'should return size -1'() {
        expect:
        new DatasetInputStream(new ByteArrayInputStream(new byte[0])).size() == -1L
    }

    def 'should throw on write'() {
        when:
        new DatasetInputStream(new ByteArrayInputStream(new byte[0])).write(ByteBuffer.allocate(1))

        then:
        thrown(UnsupportedOperationException)
    }

    def 'should throw on seek'() {
        when:
        new DatasetInputStream(new ByteArrayInputStream(new byte[0])).position(0L)

        then:
        thrown(UnsupportedOperationException)
    }

    def 'should throw on truncate'() {
        when:
        new DatasetInputStream(new ByteArrayInputStream(new byte[0])).truncate(0L)

        then:
        thrown(UnsupportedOperationException)
    }
}
