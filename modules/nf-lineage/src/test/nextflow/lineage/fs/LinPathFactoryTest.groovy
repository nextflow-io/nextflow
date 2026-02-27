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

package nextflow.lineage.fs

import java.nio.file.Files
import java.nio.file.Path

import nextflow.Global
import nextflow.Session
import spock.lang.Specification
import spock.lang.Unroll

/**
 * LID Path Factory tests.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class LinPathFactoryTest extends Specification {

    Path tmp

    def setup() {
        tmp = Files.createTempDirectory("data")
        Global.session = Mock(Session) { getConfig()>> [workflow:[lineage:[store:[location: tmp.toString()]]]] }
    }

    def cleanup() {
        Global.session = null
        tmp.deleteDir()
    }

    def 'should create lin path' () {
        given:
        def factory = new LinPathFactory()

        expect:
        factory.parseUri('foo') == null

        when:
        def p1 = factory.parseUri('lid://12345')
        then:
        p1.toUriString() == 'lid://12345'

        when:
        def p2 = factory.parseUri('lid://12345/x/y/z')
        then:
        p2.toUriString() == 'lid://12345/x/y/z'

        when:
        def p3 = factory.parseUri('lid://12345//x///y/z//')
        then:
        p3.toUriString() == 'lid://12345/x/y/z'

        when:
        factory.parseUri('lid:///12345')
        then:
        thrown(IllegalArgumentException)
    }

    @Unroll
    def 'should convert get lid uri string' () {
        given:
        def factory = new LinPathFactory()

        when:
        def lid = LinPathFactory.create(EXPECTED)
        then:
        factory.toUriString(lid) == EXPECTED

        where:
        _ | EXPECTED
        _ | 'lid://123'
        _ | 'lid://123/a/b/c'
    }
}
