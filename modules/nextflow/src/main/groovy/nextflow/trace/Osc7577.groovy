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
import java.security.MessageDigest

import groovy.transform.CompileStatic

/**
 * Encodes Program Resume Protocol (OSC 7577) reports, revision 0.5.
 *
 * <p>A report is {@code ESC ] 7577 ; pair:pair:mac=hex ESC \}. The pairs are
 * {@code pid}, {@code time}, either {@code argv} or {@code clear}, and an
 * optional {@code cwd}. {@code argv} elements and {@code cwd} are base64 of UTF-8
 * text. The {@code mac} is the SHA-256 of the body, a colon and the terminal's
 * token, and is always the last key.
 *
 * <p>Spec: https://www.superlogical.com/rex/docs/build/program-resume
 *
 * <p>The token itself is never written. Both methods return {@code null} when
 * the report breaks a rule of the spec or does not fit its size limits, so the
 * caller sends nothing instead of a report the terminal would reject.
 */
@CompileStatic
class Osc7577 {

    static final String TOKEN_ENV = 'OSC7577_TOKEN'

    /** Minimum total length every terminal must accept, in bytes, including the OSC prefix and ST */
    static final int MAX_LENGTH = 2048

    /** Minimum number of `argv` elements every terminal must accept */
    static final int MAX_ARGV = 64

    private static final String PREFIX = '\u001B]7577;'

    private static final String ST = '\u001B\\'

    /**
     * @return a report that makes the terminal resume {@code argv}, or {@code null} if it cannot be reported
     */
    static String report(String token, long pid, long time, List<String> argv, String cwd) {
        if( !validPidAndTime(pid, time) || !validArgv(argv) )
            return null
        final body = new StringBuilder(128)
        body.append('pid=').append(pid)
        body.append(':time=').append(time)
        body.append(':argv=').append(argv.collect { String it -> base64(it) }.join(','))
        if( cwd != null ) {
            if( !validPath(cwd) )
                return null
            body.append(':cwd=').append(base64(cwd))
        }
        return seal(token, body)
    }

    /**
     * @return a report that removes the record of {@code pid}
     */
    static String clear(String token, long pid, long time) {
        if( !validPidAndTime(pid, time) )
            return null
        final body = new StringBuilder(64)
        body.append('pid=').append(pid)
        body.append(':time=').append(time)
        body.append(':clear=1')
        return seal(token, body)
    }

    private static String seal(String token, StringBuilder body) {
        if( !token )
            return null
        final mac = sha256Hex("${body}:${token}".toString())
        final result = PREFIX + body + ':mac=' + mac + ST
        return result.length() <= MAX_LENGTH ? result : null
    }

    private static boolean validPidAndTime(long pid, long time) {
        return pid > 0 && time > 0
    }

    private static boolean validArgv(List<String> argv) {
        if( !argv || argv.size() > MAX_ARGV )
            return false
        for( String it : argv ) {
            if( !validText(it) )
                return false
        }
        // the program: an absolute path, or a bare name that is not an option
        final program = argv[0]
        if( !program )
            return false
        if( program.startsWith('/') )
            return true
        return !program.contains('/') && !program.startsWith('-')
    }

    private static boolean validPath(String path) {
        return validText(path) && path.startsWith('/')
    }

    /**
     * No control characters (C0, DEL, C1) and text that is valid UTF-8, which
     * rules out unpaired surrogates.
     */
    private static boolean validText(String text) {
        if( text == null )
            return false
        for( int i = 0; i < text.length(); i++ ) {
            final ch = text.charAt(i)
            if( ch < 0x20 || ch == 0x7F || (ch >= 0x80 && ch <= 0x9F) )
                return false
        }
        return StandardCharsets.UTF_8.newEncoder().canEncode(text)
    }

    private static String base64(String text) {
        return Base64.encoder.encodeToString(text.getBytes(StandardCharsets.UTF_8))
    }

    private static String sha256Hex(String text) {
        final digest = MessageDigest.getInstance('SHA-256').digest(text.getBytes(StandardCharsets.UTF_8))
        final hex = new StringBuilder(64)
        for( byte b : digest )
            hex.append(String.format('%02x', b))
        return hex.toString()
    }

}
