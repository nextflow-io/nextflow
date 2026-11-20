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

import groovy.transform.CompileStatic
import org.fusesource.jansi.internal.CLibrary

/**
 * Encodes a Program Status Protocol report (OSC 7501, revision 0.4).
 *
 * <p>The sequence is {@code ESC ] 7501 ; key=value:key=value ESC \}.
 * {@code state} is required. {@code app} is repeated on every report because
 * a report replaces its record. {@code msg} is base64 of one UTF-8 line.
 * Free text is stripped of control characters so a terminal does not discard
 * the report. See {@code docs/reference/osc-7501.mdx}.
 */
@CompileStatic
class Osc7501Emitter implements TerminalActivityEmitter {

    static final String APP = 'nextflow'

    private static final char ESC = 0x1B as char

    private final Appendable out

    Osc7501Emitter() {
        this(System.out)
    }

    Osc7501Emitter(Appendable out) {
        this.out = out
    }

    @Override
    void emit(String state, String message) {
        final text = format(state, message)
        synchronized (out) {
            out.append(text)
            if( out instanceof Flushable )
                ((Flushable) out).flush()
        }
    }

    /**
     * @return the OSC 7501 report, including the leading ESC and the ST ({@code ESC \}) terminator
     */
    static String format(String state, String message) {
        final body = new StringBuilder(80)
        body.append('state=').append(state)
        body.append(':app=').append(APP)
        final encoded = encodeMessage(message)
        if( encoded )
            body.append(':msg=').append(encoded)
        final seq = new StringBuilder(body.length() + 8)
        seq.append(ESC).append(']7501;').append(body).append(ESC).append('\\')
        return seq.toString()
    }

    /**
     * True when file descriptor 1 is a terminal.
     *
     * <p>The compile classpath's jansi is the copy inside jline 2.9. Its
     * {@code CLibrary.isatty} reports on stdout alone. {@link System#console()}
     * also requires a console on stdin, and on JDK 22–24 it can be non-null
     * without a tty, so it is only the fallback when the native check is absent.
     */
    static boolean stdoutIsTty() {
        try {
            if( CLibrary.HAVE_ISATTY )
                return CLibrary.isatty(CLibrary.STDOUT_FILENO) != 0
        }
        catch( Throwable ignored ) {
        }
        return System.console() != null
    }

    private static String encodeMessage(String message) {
        final text = sanitize(message)
        if( !text )
            return null
        return Base64.encoder.encodeToString(text.getBytes(StandardCharsets.UTF_8))
    }

    /**
     * Drop C0 and C1 controls, including newline and tab. The protocol
     * discards a report whose decoded message contains one.
     */
    private static String sanitize(String message) {
        if( !message )
            return null
        final buf = new StringBuilder(message.length())
        for( int i = 0; i < message.length(); i++ ) {
            final ch = message.charAt(i)
            if( ch >= 0x20 && ch != 0x7F && !(ch >= 0x80 && ch <= 0x9F) )
                buf.append(ch)
        }
        final text = buf.toString().trim()
        return text ? text : null
    }

}
