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

import spock.lang.Specification
import spock.lang.Unroll

/**
 * The reports and hashes below are the published examples of OSC 7577
 * revision 0.5, which use the token {@code dummy-token-replace-me}.
 */
class Osc7577EmitterTest extends Specification {

    static final String TOKEN = 'dummy-token-replace-me'

    static final long PID = 74877

    static final long TIME = 1795543200

    static final String ESC = '\u001b'

    def 'should write the neovim command with its working directory'() {
        expect:
        Osc7577Emitter.report(PID, TIME, ['nvim', '-S', '/home/user/.local/state/nvim/sessions/8f2c1d.vim'], '/home/user/src/project', TOKEN) ==
                ESC + ']7577;pid=74877:time=1795543200' +
                ':argv=bnZpbQ==,LVM=,L2hvbWUvdXNlci8ubG9jYWwvc3RhdGUvbnZpbS9zZXNzaW9ucy84ZjJjMWQudmlt' +
                ':cwd=L2hvbWUvdXNlci9zcmMvcHJvamVjdA==' +
                ':mac=cd283d150dc1e44ac3137dedbe50bcb9ac73666a4f9ee9a12f39f6adffbab247' +
                ESC + '\\'
    }

    @Unroll
    def 'should compute the published mac for #name'() {
        expect:
        Osc7577Emitter.mac(body, TOKEN) == expected

        where:
        name       | body                                                                                                                                         | expected
        'neovim'   | 'pid=74877:time=1795543200:argv=bnZpbQ==,LVM=,L2hvbWUvdXNlci8ubG9jYWwvc3RhdGUvbnZpbS9zZXNzaW9ucy84ZjJjMWQudmlt'                            | '9dbb99bfb856749a87563839d7ce9062818d78cf25909db93f309e0cd3b1909b'
        'clear'    | 'pid=74877:time=1795543200:clear=1'                                                                                                          | 'c30603369f98fef581a4098be889175b0c74f0641e49f8b25f4326cb38f97c93'
        'lazygit'  | 'pid=74877:time=1795543200:argv=bGF6eWdpdA==:cwd=L2hvbWUvdXNlci9zcmMvcHJvamVjdA=='                                                          | '35cb909805c380336efd3689a125acc774a9dec8422a1d21492ebec06ca758e1'
        'wrapper'  | 'pid=90412:time=1795543200:argv=cmVzdW1lLXdyYXA=,c3No,ZGV2Ym94'                                                                              | 'e76149d0ba54a83570ffab0583867dacdfa23928ba9d18604f0badea88b62b65'
    }

    def 'should match the published clear report'() {
        expect:
        Osc7577Emitter.clear(PID, TIME, TOKEN) ==
                ESC + ']7577;pid=74877:time=1795543200:clear=1:mac=c30603369f98fef581a4098be889175b0c74f0641e49f8b25f4326cb38f97c93' + ESC + '\\'
    }

    def 'should end the body with the mac and never write the token'() {
        when:
        def seq = Osc7577Emitter.report(PID, TIME, ['nextflow', 'run', 'main.nf'], '/work', TOKEN)

        then:
        seq.startsWith(ESC + ']7577;pid=74877:time=1795543200:argv=')
        seq.substring(seq.lastIndexOf(':mac=')) ==~ ':mac=[0-9a-f]{64}' + ESC + '\\\\'
        !seq.contains(TOKEN)
    }

    def 'should use the injected pid and clock when writing'() {
        given:
        def sink = new StringBuilder()
        def emitter = new Osc7577Emitter(sink, TOKEN, PID, { TIME })

        when:
        emitter.report(['nextflow', 'run', 'main.nf'], '/work')
        emitter.clear()

        then:
        sink.toString() ==
                Osc7577Emitter.report(PID, TIME, ['nextflow', 'run', 'main.nf'], '/work', TOKEN) +
                Osc7577Emitter.clear(PID, TIME, TOKEN)
    }

    def 'should encode non ascii arguments as utf-8 base64'() {
        when:
        def seq = Osc7577Emitter.report(PID, TIME, ['nextflow', 'run', 'caf\u00e9.nf'], '/work', TOKEN)

        then:
        seq.contains(Base64.encoder.encodeToString('caf\u00e9.nf'.getBytes('UTF-8')))
    }

    @Unroll
    def 'should accept the token #token'() {
        expect:
        Osc7577Emitter.isValidToken(token) == valid

        where:
        token                                  | valid
        'dummy-token-replace-me'               | true
        'A' * 64                               | true
        'A' * 21                               | false
        'A' * 65                               | false
        'dummy-token-replace-m!'               | false
        'dummy token replace me'               | false
        null                                   | false
        ''                                     | false
    }

    def 'should refuse an invalid token'() {
        when:
        new Osc7577Emitter(new StringBuilder(), 'short', PID, { TIME })

        then:
        thrown(IllegalArgumentException)
    }

    @Unroll
    def 'should refuse a command the protocol cannot carry: #name'() {
        when:
        Osc7577Emitter.report(PID, TIME, argv, cwd, TOKEN)

        then:
        thrown(IllegalArgumentException)

        where:
        name                          | argv                                 | cwd
        'empty command'               | []                                   | '/work'
        'relative program path'       | ['./run.sh']                         | '/work'
        'program starting with dash'  | ['-nextflow']                        | '/work'
        'empty program'               | ['']                                 | '/work'
        'newline in an argument'      | ['nextflow', 'a\nb']                 | '/work'
        'escape in an argument'       | ['nextflow', 'a\u001bb']             | '/work'
        'c1 control in an argument'   | ['nextflow', 'a\u009bb']             | '/work'
        'delete in an argument'       | ['nextflow', 'a\u007fb']             | '/work'
        'lone surrogate'              | ['nextflow', 'a\ud800b']             | '/work'
        'relative cwd'                | ['nextflow']                         | 'work'
        'missing cwd'                 | ['nextflow']                         | null
        'too many arguments'          | ['nextflow'] + (1..64).collect { 'x' } | '/work'
        'longer than 2048 bytes'      | ['nextflow', 'x' * 2000]             | '/work'
    }

    def 'should allow an absolute program path and empty arguments after the first'() {
        when:
        def seq = Osc7577Emitter.report(PID, TIME, ['/usr/local/bin/nextflow', '', 'run'], '/work', TOKEN)

        then:
        seq.contains('argv=L3Vzci9sb2NhbC9iaW4vbmV4dGZsb3c=,,cnVu:')
    }

    def 'should allow exactly 64 arguments'() {
        when:
        Osc7577Emitter.report(PID, TIME, ['nextflow'] + (1..63).collect { 'x' }, '/work', TOKEN)

        then:
        noExceptionThrown()
    }

}
