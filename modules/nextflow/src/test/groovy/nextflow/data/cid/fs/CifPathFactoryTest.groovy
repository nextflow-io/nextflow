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

package nextflow.data.cid.fs

import java.nio.file.Path

import nextflow.Global
import nextflow.Session
import spock.lang.Specification
import spock.lang.Unroll

/**
 * CID Path Factory tests.
 *
 * @author Jorge Ejarque <jorge.ejarque@seqera.io>
 */
class CifPathFactoryTest extends Specification {

    def setup() {
        Global.session = Mock(Session) { getConfig()>> [workflow:[data:[store:[location: '/some/data']]]] }
    }

    def cleanup() {
        Global.session = null
    }

    def 'should create cid path' () {
        given:
        def factory = new CidPathFactory()

        expect:
        factory.parseUri('foo') == null
        
        when:
        def p1 = factory.parseUri('cid://12345')
        then:
        p1.getTargetPath()  == Path.of('/some/data/.meta/12345')
        p1.toUriString() == 'cid://12345'

        when:
        def p2 = factory.parseUri('cid://12345/x/y/z')
        then:
        p2.getTargetPath() == Path.of('/some/data/.meta/12345/x/y/z')
        p2.toUriString() == 'cid://12345/x/y/z'

        when:
        def p3 = factory.parseUri('cid://12345//x///y/z//')
        then:
        p3.getTargetPath() == Path.of('/some/data/.meta/12345/x/y/z')
        p2.toUriString() == 'cid://12345/x/y/z'

        when:
        factory.parseUri('cid:///12345')
        then:
        thrown(IllegalArgumentException)
    }

    @Unroll
    def 'should convert get cid uri string' () {
        given:
        def factory = new CidPathFactory()

        when:
        def cid = CidPathFactory.create(EXPECTED)
        then:
        factory.toUriString(cid) == EXPECTED
        
        where:
        _ | EXPECTED
        _ | 'cid://123'
        _ | 'cid://123/a/b/c'
    }
}
