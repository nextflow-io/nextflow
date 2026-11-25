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

import java.nio.CharBuffer
import java.nio.charset.CharacterCodingException
import java.nio.charset.CodingErrorAction
import java.nio.charset.StandardCharsets
import java.security.MessageDigest
import java.util.function.LongSupplier
import java.util.regex.Pattern

import groovy.transform.CompileStatic

/**
 * Encodes Program Resume Protocol reports (OSC 7577, revision 0.5).
 *
 * <p>The sequence is {@code ESC ] 7577 ; pid=..:time=..:argv=..[:cwd=..]:mac=.. ESC \}.
 * The body is {@code key=value} pairs joined by {@code :}, and {@code mac} is
 * always last. {@code mac} is the SHA-256 of the body before {@code :mac=},
 * then {@code :}, then the token the terminal put in {@code OSC7577_TOKEN}.
 * The token itself is never written.
 *
 * <p>A report that the protocol cannot carry (control characters, a relative
 * program path, more than 64 arguments, more than 2048 bytes) is refused with
 * an {@link IllegalArgumentException} instead of being sent in a form the
 * terminal would reject.
 *
 * @see <a href="https://www.superlogical.com/rex/docs/build/program-resume">Program Resume Protocol</a>
 */
@CompileStatic
class Osc7577Emitter {

    /** Minimum total sequence length every terminal must accept, in bytes */
    static final int MAX_LENGTH = 2048

    /** Minimum number of {@code argv} elements every terminal must accept */
    static final int MAX_ARGV = 64

    private static final char ESC = 0x1B as char

    private static final Pattern TOKEN = ~/[A-Za-z0-9_-]{22,64}/

    private final Appendable out

    private final String token

    private final LongSupplier epochSeconds

    private final long pid

    Osc7577Emitter(Appendable out, String token, long pid, LongSupplier epochSeconds) {
        if( !isValidToken(token) )
            throw new IllegalArgumentException("Not a valid OSC 7577 token")
        this.out = out
        this.token = token
        this.pid = pid
        this.epochSeconds = epochSeconds
    }

    /**
     * Ask the terminal to bring {@code argv} back after a restart.
     *
     * @param argv the resume command, one element per argument
     * @param cwd absolute working directory to resume in
     */
    void report(List<String> argv, String cwd) {
        write(report(pid, epochSeconds.getAsLong(), argv, cwd, token))
    }

    /** Ask the terminal to forget the resume command of this process. */
    void clear() {
        write(clear(pid, epochSeconds.getAsLong(), token))
    }

    private void write(String text) {
        synchronized (out) {
            out.append(text)
            if( out instanceof Flushable )
                ((Flushable) out).flush()
        }
    }

    static boolean isValidToken(String token) {
        return token != null && TOKEN.matcher(token).matches()
    }

    /**
     * @return the resume report, including the leading ESC and the ST ({@code ESC \}) terminator
     * @throws IllegalArgumentException when the protocol cannot carry the command
     */
    static String report(long pid, long time, List<String> argv, String cwd, String token) {
        validateArgv(argv)
        validateAbsolute(cwd, 'cwd')
        final body = new StringBuilder(256)
        body.append('pid=').append(pid)
        body.append(':time=').append(time)
        body.append(':argv=').append(argv.collect { String it -> encode(it) }.join(','))
        body.append(':cwd=').append(encode(cwd))
        return seal(body.toString(), token)
    }

    /**
     * @return the report that removes the record of {@code pid}, including ESC and ST
     */
    static String clear(long pid, long time, String token) {
        return seal("pid=${pid}:time=${time}:clear=1".toString(), token)
    }

    private static String seal(String body, String token) {
        final seq = new StringBuilder(body.length() + 96)
        seq.append(ESC).append(']7577;').append(body)
        seq.append(':mac=').append(mac(body, token))
        seq.append(ESC).append('\\')
        if( seq.length() > MAX_LENGTH )
            throw new IllegalArgumentException("Resume report is ${seq.length()} bytes, the protocol guarantees ${MAX_LENGTH}")
        return seq.toString()
    }

    /** SHA-256 over {@code body}, {@code :}, and the token, as 64 lowercase hex characters */
    static String mac(String body, String token) {
        final digest = MessageDigest.getInstance('SHA-256')
        final hash = digest.digest("${body}:${token}".toString().getBytes(StandardCharsets.UTF_8))
        final hex = new StringBuilder(64)
        for( byte b : hash )
            hex.append(String.format('%02x', b & 0xFF))
        return hex.toString()
    }

    private static void validateArgv(List<String> argv) {
        if( !argv )
            throw new IllegalArgumentException('Resume command is empty')
        if( argv.size() > MAX_ARGV )
            throw new IllegalArgumentException("Resume command has ${argv.size()} arguments, the protocol guarantees ${MAX_ARGV}")
        final program = argv[0]
        if( !program )
            throw new IllegalArgumentException('Resume command has an empty program name')
        if( program.startsWith('-') )
            throw new IllegalArgumentException("Program name cannot start with '-': ${program}")
        if( program.contains('/') )
            validateAbsolute(program, 'program path')
        for( String it : argv )
            checkText(it, 'argument')
    }

    private static void validateAbsolute(String path, String what) {
        if( !path || !path.startsWith('/') )
            throw new IllegalArgumentException("Resume ${what} must be absolute: ${path}")
        checkText(path, what)
    }

    private static void checkText(String text, String what) {
        for( int i = 0; i < text.length(); i++ ) {
            final ch = text.charAt(i)
            if( ch < 0x20 || ch == 0x7F || (ch >= 0x80 && ch <= 0x9F) )
                throw new IllegalArgumentException("Resume ${what} contains a control character")
        }
    }

    private static String encode(String text) {
        try {
            final bytes = StandardCharsets.UTF_8.newEncoder()
                    .onMalformedInput(CodingErrorAction.REPORT)
                    .onUnmappableCharacter(CodingErrorAction.REPORT)
                    .encode(CharBuffer.wrap(text))
            final copy = new byte[bytes.remaining()]
            bytes.get(copy)
            return Base64.encoder.encodeToString(copy)
        }
        catch( CharacterCodingException e ) {
            throw new IllegalArgumentException('Resume text is not valid UTF-8', e)
        }
    }

}
