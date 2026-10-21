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

package nextflow.util

import spock.lang.Specification
import spock.lang.Timeout
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
        null                                    | null
        [foo:'Hello']                           | [foo:'Hello']
        [foo: [bar: 'World']]                   | [foo: [bar: 'World']]
        [foo: [password:'hola', token:'hi']]    | [foo: [password:'****', token:'****']]
        [foo: [password:'1234567890']]          | [foo: [password:'****']]
        [foo: [password:null, token:'']]        | [foo: [password:'(null)', token:'(empty)']]
        [foo: [customPassword:'hola']]          | [foo: [customPassword:'****']]
        [foo: [towerLicense:'hola']]            | [foo: [towerLicense:'****']]
        [url: 'redis://host:port']              | [url: 'redis://host:port']
        [url: 'redis://secret@host:port']       | [url: 'redis://****@host:port']
        [url: 'ftp://secretlong@host:port/x/y'] | [url: 'ftp://sec****@host:port/x/y']
        [providers:[github:[auth:'12345']]]     | [providers:[github:[auth:'****']]]
        [AZURE_STORAGE_ACCOUNT_KEY:'1234567890'] | [AZURE_STORAGE_ACCOUNT_KEY:'****']
        [sortKey:'name', region:'eu']           | [sortKey:'name', region:'eu']
        [basicAuth:'x', authorization:'y', tokens:'z'] | [basicAuth:'****', authorization:'****', tokens:'****']
        [azure:[storage:[accountKey:'xyz']]]    | [azure:[storage:[accountKey:'****']]]
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
        'providers.github.auth=\'user:ghp_a-b\''        | 'providers.github.auth=\'********\''
        'docker run -e "AWS_SECRET_ACCESS_KEY=a/b+c" -e "FUSION_LICENSE_TOKEN=x.y-z" -e "FOO=bar"' | 'docker run -e "AWS_SECRET_ACCESS_KEY=********" -e "FUSION_LICENSE_TOKEN=********" -e "FOO=bar"'
        '[AZURE_STORAGE_ACCOUNT_KEY:abc==, NAME:foo]'   | '[AZURE_STORAGE_ACCOUNT_KEY:********, NAME:foo]'
        'token: "********"'                             | 'token: "********"'
        'password=*abc123'                              | 'password=********'
        'password=ab*cdef'                              | 'password=********'
        'password=\'a,b\''                              | 'password=\'********\''
        'password="my secret pass"'                     | 'password="********"'
        'token    = \'ghp_xxx\''                        | 'token    = \'********\''
        'password:\t\t"abc"'                             | 'password:\t\t"********"'
        'SINGULARITYENV_AWS_SESSION_TOKEN="a b"'        | 'SINGULARITYENV_AWS_SESSION_TOKEN="********"'
        'password:\nabc'                                 | 'password:\n********'
        'password\n=abc'                                 | 'password\n=********'
        'password="abc'                                  | 'password="********'
        'password=\'abc'                                 | 'password=\'********'
        'Authorization: Bearer abc.def'                 | 'Authorization: ********'
        'authorization=abc'                             | 'authorization=********'
        'oauth=abc x_auth=def'                          | 'oauth=******** x_auth=********'
        'X-Api-Key: abc'                                | 'X-Api-Key: ********'
        'api-key=abc'                                   | 'api-key=********'
        'access_tokens=abc'                             | 'access_tokens=********'
        'password=abc,def'                              | 'password=********'
        'password: \'a,b'                               | 'password: \'********'
        '[A_TOKEN:x, NAME:foo]'                         | '[A_TOKEN:********, NAME:foo]'
        // not secrets
        'monkey: banana, Author: Paolo'                 | 'monkey: banana, Author: Paolo'
        'sortKey=name tokens: 5'                        | 'sortKey=name tokens: 5'
        '-v /home/u/.ssh/authorized_keys:/root/x'       | '-v /home/u/.ssh/authorized_keys:/root/x'
        'quay.io/x/hotkey:1.0'                          | 'quay.io/x/hotkey:1.0'
        '-e "FUSION_AWS_SSEKMS_KEY_ID=arn:aws:kms:x"'   | '-e "FUSION_AWS_SSEKMS_KEY_ID=arn:aws:kms:x"'
        'cost: $1 \\ x'                                 | 'cost: $1 \\ x'
    }

    @Unroll
    def 'should check sensitive name' () {
        expect:
        StringUtils.isSensitiveName(NAME) == EXPECTED

        where:
        NAME                        | EXPECTED
        'password'                  | true
        'DB_PASSWORD'               | true
        'passwd'                    | true
        'accessToken'               | true
        'tokens'                    | true
        'towerLicense'              | true
        'auth'                      | true
        'basicAuth'                 | true
        'credentials'               | true
        'secretKey'                 | true
        'accessKey'                 | true
        'accountKey'                | true
        'AZURE_STORAGE_ACCOUNT_KEY' | true
        'apiKey'                    | true
        'OPENAI_API_KEY'            | true
        'privateKey'                | true
        'AWS_ACCESS_KEY_ID'         | true
        'AWS_REGION'                | true
        'API_KEY'                   | true
        and:
        'user'                      | false
        'region'                    | false
        'aws'                       | false
        'sortKey'                   | false
        'monkey'                    | false
        'FUSION_AWS_SSEKMS_KEY_ID'  | false
        ''                          | false
        null                        | false
    }

    @Timeout(5)
    def 'should strip secrets of long strings in linear time' () {
        expect:
        [ 'a' * 100_000, 'key' * 30_000, 'password' * 10_000, 'x=' * 50_000, 'a-' * 50_000 ].every { StringUtils.stripSecrets(it) == it }
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
