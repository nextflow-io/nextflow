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

package nextflow.util

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class StringUtilsTest extends Specification {

    @Unroll
    def 'should get url protocol' () {
        expect:
        StringUtils.getUrlProtocol(STR)  == EXPECTED
        where:
        EXPECTED    | STR
        'ftp'       | 'ftp://abc.com'
        's3'        | 's3://bucket/abc'
        null        | 'gitzabc:xyz'
        null        | '/a/bc/'
    }

    @Unroll
    def 'should find base url #STR' () {
        expect:
        StringUtils.baseUrl(STR) == BASE

        where:
        BASE                | STR
        null                | null
        'http://foo.com'    | 'http://foo.com'
        'http://foo.com'    | 'http://foo.com/abc'
        'http://foo.com'    | 'http://foo.com/abc/mskd0fs =ds0f'
        and:
        'https://foo.com'    | 'https://foo.com'
        'https://foo.com'    | 'https://foo.com/abc'
        'https://foo.com'    | 'https://foo.com/abc/mskd0fs =ds0f'
        and:
        'https://foo.com'    | 'HTTPS://FOO.COM'
        'https://foo.com'    | 'HTTPS://FOO.COM/ABC'
        'https://foo.com'    | 'HTTPS://FOO.COM/ABC/MSKD0FS =DS0F'
        and:
        'https://foo.com:80' | 'https://foo.com:80'
        'https://foo.com:80' | 'https://foo.com:80/abc'
        'https://foo.com:80' | 'https://foo.com:80/abc/mskd0fs =ds0f'
        and:
        'ftp://foo.com:80'   | 'ftp://foo.com:80'
        'ftp://foo.com:80'   | 'ftp://foo.com:80/abc'
        's3://foo.com:80'    | 's3://foo.com:80/abc'
        and:
        'dx://project-123'   | 'dx://project-123'
        'dx://project-123:'  | 'dx://project-123:'
        'dx://project-123:'  | 'dx://project-123:/abc'
        and:
        null                 | 'blah'
        null                 | 'http:/xyz.com'
        null                 | '1234://xyz'
        null                 | '1234://xyz.com/abc'
    }

    @Unroll
    def 'should strip passwords' () {
        expect:
        StringUtils.stripSecrets(SECRET) == EXPECTED

        where:
        SECRET                                  | EXPECTED
        [foo:'Hello']                           | [foo:'Hello']
        [foo: [bar: 'World']]                   | [foo: [bar: 'World']]
        [foo: [password:'hola', token:'hi']]    | [foo: [password:'****', token:'****']]
        [foo: [password:'1234567890']]          | [foo: [password:'123****']]
        [foo: [customPassword:'hola']]          | [foo: [customPassword:'****']]
        [foo: [towerLicense:'hola']]            | [foo: [towerLicense:'****']]
        [url: 'redis://host:port']              | [url: 'redis://host:port']
        [url: 'redis://secret@host:port']       | [url: 'redis://****@host:port']
        [url: 'ftp://secretlong@host:port/x/y'] | [url: 'ftp://sec****@host:port/x/y']
        [providers:[github:[auth:'12345']]]     | [providers:[github:[auth:'****']]]
    }

    @Unroll
    def 'should strip sensitive strings' () {
        expect:
        StringUtils.stripSecrets(SECRET) == EXPECTED

        where:
        SECRET                                          | EXPECTED
        'Hi\n here is the "password" : "1234"'          | 'Hi\n here is the "password" : "********"'
        'Hi\n here is the password : "1"'               | 'Hi\n here is the password : "********"'
        'Hi\n here is the password : \'1\''             | 'Hi\n here is the password : \'********\''
        'Hi\n "password" :"1" \n "token": "123"'        | 'Hi\n "password" :"********" \n "token": "********"'
        'Hi\n "password" :\'1\' \n "token": "123"'      | 'Hi\n "password" :\'********\' \n "token": "********"'
        'Hi\n \'password\' :\'1\' \n \'token\': \'123\''| 'Hi\n \'password\' :\'********\' \n \'token\': \'********\''
        'Hi\n password :"1"\nsecret: "345"'             | 'Hi\n password :"********"\nsecret: "********"'
        'secret="abc" password:"1" more text'           | 'secret="********" password:"********" more text'
    }

    @Unroll
    def 'should strip secret' () {
        expect:
        StringUtils.redact(SECRET) == EXPECTED

        where:
        SECRET          | EXPECTED
        'hi'            | '****'
        'Hello'         | '****'
        'World'         | '****'
        '1234567890'    | '123****'
        'hola'          | '****'
        null            | '(null)'
        ''              | '(empty)'
    }

    @Unroll
    def 'should strip url password' () {
        expect:
        StringUtils.redactUrlPassword(SECRET) == EXPECTED

        where:
        SECRET                  | EXPECTED
        'hi'                    | 'hi'
        'http://foo/bar'        | 'http://foo/bar'
        'http://secret@foo/bar' | 'http://****@foo/bar'
    }

    def 'should check ipv6' () {
        expect:
        StringUtils.isIpV6String('2001:db8:42:1:2:3:1:443')
    }

    @Unroll
    def 'should check is ipv6'  () {
        expect:
        StringUtils.isIpV6String(ADDR) == EXPECTED
        where:
        ADDR                    | EXPECTED
        null                    | false
        'foo.com'               | false
        'foo.com'               | false
        '127.0.0.0'             | false
        '2001:db8:42:1::1:443'  | true
    }

    @Unroll
    def 'should format host name'  () {
        expect:
        StringUtils.formatHostName(HOST, PORT) == EXPECTED
        where:
        HOST                    | PORT      | EXPECTED
        'foo.com'               | null      | 'foo.com'
        'foo.com'               | '80'      | 'foo.com:80'
        '127.0.0.0'             | '8000'    | '127.0.0.0:8000'
        '2001:db8:42:1::1:443'  | null      | '2001:db8:42:1::1:443'
        '2001:db8:42:1::1:443'  | '8000'    | '[2001:db8:42:1::1:443]:8000'
    }
}
