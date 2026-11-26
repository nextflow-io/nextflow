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
 * The expected reports are the examples published in the OSC 7577 specification
 * (revision 0.5), which use the token {@code dummy-token-replace-me}.
 */
class Osc7577Test extends Specification {

    static final String TOKEN = 'dummy-token-replace-me'

    static final String OSC = '\u001B]7577;'

    static final String ST = '\u001B\\'

    def 'should match the spec example for neovim' () {
        expect:
        Osc7577.report(TOKEN, 74877, 1795543200, ['nvim', '-S', '/home/user/.local/state/nvim/sessions/8f2c1d.vim'], null) ==
            OSC + 'pid=74877:time=1795543200:argv=bnZpbQ==,LVM=,L2hvbWUvdXNlci8ubG9jYWwvc3RhdGUvbnZpbS9zZXNzaW9ucy84ZjJjMWQudmlt:mac=9dbb99bfb856749a87563839d7ce9062818d78cf25909db93f309e0cd3b1909b' + ST
    }

    def 'should match the spec example for clear' () {
        expect:
        Osc7577.clear(TOKEN, 74877, 1795543200) ==
            OSC + 'pid=74877:time=1795543200:clear=1:mac=c30603369f98fef581a4098be889175b0c74f0641e49f8b25f4326cb38f97c93' + ST
    }

    def 'should match the spec example with a working directory' () {
        expect:
        Osc7577.report(TOKEN, 74877, 1795543200, ['lazygit'], '/home/user/src/project') ==
            OSC + 'pid=74877:time=1795543200:argv=bGF6eWdpdA==:cwd=L2hvbWUvdXNlci9zcmMvcHJvamVjdA==:mac=35cb909805c380336efd3689a125acc774a9dec8422a1d21492ebec06ca758e1' + ST
    }

    def 'should match the spec example for the wrapper' () {
        expect:
        Osc7577.report(TOKEN, 90412, 1795543200, ['resume-wrap', 'ssh', 'devbox'], null) ==
            OSC + 'pid=90412:time=1795543200:argv=cmVzdW1lLXdyYXA=,c3No,ZGV2Ym94:mac=e76149d0ba54a83570ffab0583867dacdfa23928ba9d18604f0badea88b62b65' + ST
    }

    def 'should never write the token' () {
        expect:
        !Osc7577.report(TOKEN, 1, 1, ['nextflow'], '/tmp').contains(TOKEN)
        !Osc7577.clear(TOKEN, 1, 1).contains(TOKEN)
    }

    def 'should keep arguments with shell syntax as plain text' () {
        when:
        final text = Osc7577.report(TOKEN, 1, 1, ['nextflow', '$(reboot)', 'a; b', ''], null)
        final argv = (text =~ /argv=([^:]*):/)[0][1].toString().split(',', -1)
        then:
        argv.collect { new String(Base64.decoder.decode(it), 'UTF-8') } == ['nextflow', '$(reboot)', 'a; b', '']
    }

    @Unroll
    def 'should not report argv #argv' () {
        expect:
        Osc7577.report(TOKEN, 1, 1, argv, null) == null

        where:
        argv                                  | _
        null                                  | _
        []                                    | _
        ['']                                  | _
        ['./run.sh']                          | _
        ['bin/nextflow']                      | _
        ['-nextflow']                         | _
        ['nextflow', 'a\nb']                  | _
        ['nextflow', 'a\u0000b']              | _
        ['nextflow', 'a\u007Fb']              | _
        ['nextflow', 'a\u0085b']              | _
        ['nextflow', null]                    | _
        ['nextflow', '\uD800']                | _
        (['nextflow'] + (1..64).collect { 'x' }) | _
    }

    def 'should accept the most argv elements the spec guarantees' () {
        expect:
        Osc7577.report(TOKEN, 1, 1, ['nextflow'] + (1..63).collect { 'x' }, null) != null
    }

    def 'should not report a body that exceeds the total length' () {
        when:
        final fits = Osc7577.report(TOKEN, 1, 1, ['nextflow', 'x' * 700], null)
        final tooLong = Osc7577.report(TOKEN, 1, 1, ['nextflow', 'x' * 1600], null)
        then:
        fits != null
        fits.length() <= Osc7577.MAX_LENGTH
        tooLong == null
    }

    @Unroll
    def 'should not report cwd #cwd' () {
        expect:
        Osc7577.report(TOKEN, 1, 1, ['nextflow'], cwd) == null

        where:
        cwd         | _
        'project'   | _
        './project' | _
        '/tmp/a\nb' | _
    }

    @Unroll
    def 'should not report pid #pid and time #time' () {
        expect:
        Osc7577.report(TOKEN, pid, time, ['nextflow'], null) == null
        Osc7577.clear(TOKEN, pid, time) == null

        where:
        pid | time
        0   | 1
        -1  | 1
        1   | 0
        1   | -5
    }

    def 'should not report without a token' () {
        expect:
        Osc7577.report(token, 1, 1, ['nextflow'], null) == null
        Osc7577.clear(token, 1, 1) == null

        where:
        token << [null, '']
    }

    def 'should end with the mac and the string terminator' () {
        when:
        final text = Osc7577.report(TOKEN, 1, 1, ['nextflow'], '/tmp')
        then:
        text.startsWith(OSC)
        text.endsWith(ST)
        text =~ /:mac=[0-9a-f]{64}\u001B\\$/
    }

}
