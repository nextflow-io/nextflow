/*
 * Copyright 2013-2024, Seqera Labs
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

package io.seqera.wave.plugin.util

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class GnuCliOptsTest extends Specification {

    def 'should parse cli' () {

        when:
        def cli = GnuCliOpts.parse(['-i', '-v', 'x=y', 'nextflow', '--', 'this', '--that'])
        then:
        cli.options == ['-i':'','-v':'x=y']
        cli.container == 'nextflow'
        cli.args == ['this','--that']

        when:
        cli = GnuCliOpts.parse(['-v', 'x=y', '-w', '$PWD', 'nextflow', '-it', '--', 'this', '--that'])
        then:
        cli.options == ['-v':'x=y', '-w':'$PWD', '-it':'']
        cli.container == 'nextflow'
        cli.args == ['this','--that']
    }

}
