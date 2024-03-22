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

package nextflow.config.scope.nextflow


import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PublishDirOptsTest extends Specification {

    def 'should validate equals and hash code' () {
        when:
        def p1 = new PublishDirOpts(enabled: true, mode:'foo')
        def p2 = new PublishDirOpts(enabled: true, mode:'foo')
        def p3 = new PublishDirOpts(enabled: false, mode:'bar')

        then:
        p1 == p2 
        p1 != p3
        and:
        p1.hashCode() == p2.hashCode()
        p1.hashCode() != p3.hashCode()
    }

    def 'should create empty opts' () {
        given:
        def opts = new PublishDirOpts(Map.of())
        expect:
        opts.mode == null
        opts.enabled == null
        opts.failOnError == null
        opts.pattern == null
        opts.contentType == null
        opts.overwrite == null
        opts.storageClass == null
        opts.tags == null
    }

    @Unroll
    def 'should populate publishdir obj' () {
        expect:
        new PublishDirOpts(OPTS)."$PROPERTY" == EXPECTED
        where:
        OPTS                    | PROPERTY  | EXPECTED
        [mode:'foo']            | 'mode'    | 'foo'
        [mode:'bar']            | 'mode'    | 'bar'
        and:
        [:]                     | 'enabled' | null
        [enabled:false]         | 'enabled' | false
        [enabled:'false']       | 'enabled' | false
        [enabled:true]          | 'enabled' | true
        [enabled:'true']        | 'enabled' | true
        and:
        [:]                     | 'overwrite' | null
        [overwrite:false]       | 'overwrite' | false
        [overwrite:'false']     | 'overwrite' | false
        [overwrite:true]        | 'overwrite' | true
        [overwrite:'true']      | 'overwrite' | true
        and:
        [:]                     | 'failOnError' | null
        [failOnError:false]     | 'failOnError' | false
        [failOnError:'false']   | 'failOnError' | false
        [failOnError:true]      | 'failOnError' | true
        [failOnError:'true']    | 'failOnError' | true
        and:
        [pattern: '*.txt']      | 'pattern'     | '*.txt'
        and:
        [tags: [FOO:'one',BAR:'two']]       | 'tags'    | [FOO:'one',BAR:'two']
    }
}
