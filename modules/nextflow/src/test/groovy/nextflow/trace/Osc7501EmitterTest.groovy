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
package nextflow.trace

import java.nio.charset.StandardCharsets

import spock.lang.Specification
import spock.lang.Unroll

/**
 * Wire format for OSC 7501 revision 0.4.
 */
class Osc7501EmitterTest extends Specification {

    static final char ESC = 0x1B as char

    def 'should encode each program status the protocol defines for a run'() {
        expect:
        Osc7501Emitter.format(STATE, MESSAGE) == EXPECTED
        and:
        def seq = EXPECTED
        seq.charAt(0) == ESC
        seq.startsWith("${ESC}]7501;")
        seq.charAt(seq.length() - 2) == ESC
        seq.charAt(seq.length() - 1) == '\\' as char
        !seq.contains(MESSAGE)

        where:
        STATE      | MESSAGE                | EXPECTED
        'idle'     | 'Waiting for tasks'    | "${ESC}]7501;state=idle:app=nextflow:msg=V2FpdGluZyBmb3IgdGFza3M=${ESC}\\"
        'working'  | 'Running tasks'        | "${ESC}]7501;state=working:app=nextflow:msg=UnVubmluZyB0YXNrcw==${ESC}\\"
        'done'     | 'Workflow completed'   | "${ESC}]7501;state=done:app=nextflow:msg=V29ya2Zsb3cgY29tcGxldGVk${ESC}\\"
        'error'    | 'Workflow failed'      | "${ESC}]7501;state=error:app=nextflow:msg=V29ya2Zsb3cgZmFpbGVk${ESC}\\"
        'idle'     | 'Run cancelled'        | "${ESC}]7501;state=idle:app=nextflow:msg=UnVuIGNhbmNlbGxlZA==${ESC}\\"
    }

    def 'should omit an empty message and still name the program'() {
        expect:
        Osc7501Emitter.format('idle', null) == "${ESC}]7501;state=idle:app=nextflow${ESC}\\"
        Osc7501Emitter.format('idle', '') == "${ESC}]7501;state=idle:app=nextflow${ESC}\\"
        Osc7501Emitter.format('idle', " \n\t") == "${ESC}]7501;state=idle:app=nextflow${ESC}\\"
    }

    def 'should keep colons and controls out of the pair list'() {
        given:
        def seq = Osc7501Emitter.format('working', "a:b;c\nnext\u0007")
        def body = seq.substring(7, seq.length() - 2)
        def encoded = body.substring(body.indexOf('msg=') + 4)
        def decoded = new String(Base64.decoder.decode(encoded), StandardCharsets.UTF_8)

        expect:
        body == "state=working:app=nextflow:msg=${encoded}"
        body.count(':') == 2
        !body.contains(';')
        !body.contains('a:b')
        decoded == 'a:b;cnext'
    }

    def 'should write the report and flush'() {
        given:
        def sink = new FlushTrackingBuffer()
        def emitter = new Osc7501Emitter(sink)

        when:
        emitter.emit('done', 'Workflow completed')

        then:
        sink.toString() == Osc7501Emitter.format('done', 'Workflow completed')
        sink.flushed
    }

    @Unroll
    def 'should decode the published message for #state'() {
        given:
        def seq = Osc7501Emitter.format(state, message)
        def encoded = seq.substring(seq.indexOf('msg=') + 4, seq.length() - 2)

        expect:
        new String(Base64.decoder.decode(encoded), StandardCharsets.UTF_8) == message

        where:
        state      | message
        'idle'     | 'Waiting for tasks'
        'working'  | 'Running tasks'
        'done'     | 'Workflow completed'
        'error'    | 'Workflow failed'
    }

    private static class FlushTrackingBuffer implements Appendable, Flushable {
        private final StringBuilder buf = new StringBuilder()
        boolean flushed

        @Override
        Appendable append(CharSequence csq) {
            buf.append(csq)
            return this
        }

        @Override
        Appendable append(CharSequence csq, int start, int end) {
            buf.append(csq, start, end)
            return this
        }

        @Override
        Appendable append(char c) {
            buf.append(c)
            return this
        }

        @Override
        void flush() { flushed = true }

        @Override
        String toString() { buf.toString() }
    }

}
