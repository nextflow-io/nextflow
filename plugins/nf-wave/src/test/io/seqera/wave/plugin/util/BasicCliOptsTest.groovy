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
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BasicCliOptsTest extends Specification {

    @Unroll
    def 'should parse options' () {
        when:
        def cli = BasicCliOpts.parse(CLI?.tokenize(' '))
        then:
        cli.options == OPTS
        cli.args == ARGS
        
        where:
        CLI                                 | OPTS                                  | ARGS
        null                                | [:]                                   | []
        ''                                  | [:]                                   | []
        'alpha=1 beta=2 foo delta=3 bar'    | [alpha:'1', beta:'2', delta:'3']      | ['foo','bar']
        'alpha= foo'                        | [alpha:'']                            | ['foo']
    }
}
