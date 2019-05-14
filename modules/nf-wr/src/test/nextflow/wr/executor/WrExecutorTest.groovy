/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.executor

import nextflow.Session
import spock.lang.Specification
import java.io.BufferedWriter
import java.io.FileWriter
import nextflow.wr.client.WrRestApi

/**
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * Based on TesExecutorTest by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WrExecutorTest extends Specification {

    def 'should get endpoint' () {
        given:
        def session = Mock(Session)
        def exec = new WrExecutor(session: session)

        when:
        def result = exec.getEndPoint()
        then:
        session.getConfigAttribute('executor.wr.endpoint',_) >> 'https://localhost:11303'
        result == 'https://localhost:11303'
    }

    def 'should resolve endpoint from env'() {
        given:
        def ENV = [NXF_EXECUTOR_WR_ENDPOINT: 'https://back.end.com:11303' ]
        def session = Spy(Session)
        def exec = Spy(WrExecutor)

        when:
        def result = exec.getEndPoint()
        then:
        1 * exec.getSession() >> session
        1 * session.getSystemEnv() >> ENV
        result == 'https://back.end.com:11303'
    }

    def 'should resolve endpoint based on uid'() {
        given:
        def session = Spy(Session)
        def exec = Spy(WrExecutor)

        when:
        def result = exec.getEndPoint()
        then:
        1 * exec.getSession() >> session
        1 * exec.getUID() >> 4
        result == 'https://localhost:1038'
    }

    def 'should resolve endpoint based on deployment'() {
        given:
        def session = Spy(Session)
        def exec = Spy(WrExecutor)

        when:
        exec.resolveDeployment()
        def result = exec.getEndPoint()
        then:
        exec.getSession() >> session
        1 * session.getConfigAttribute('executor.wr.deployment',_) >> 'development'
        1 * exec.getUID() >> 4
        result == 'https://localhost:1040'
    }

    def 'should get cacert path' () {
        given:
        def session = Mock(Session)
        def exec = new WrExecutor(session: session)

        when:
        def result = exec.getCacertPath()
        then:
        session.getConfigAttribute('executor.wr.cacertpath',_) >> '/certs/ca.pem'
        result == '/certs/ca.pem'
    }

    def 'should resolve cacert path based on deployment'() {
        given:
        def session = Spy(Session)
        def exec = Spy(WrExecutor)

        when:
        exec.resolveDeployment()
        exec.resolveManagerDir()
        def result = exec.getCacertPath()
        then:
        exec.getSession() >> session
        1 * session.getConfigAttribute('executor.wr.deployment',_) >> 'production'
        1 * exec.getHomeDir() >> '/home/user'
        result == '/home/user/.wr_production/ca.pem'
    }

    def 'should get token and client'() {
        given:
        def session = Mock(Session)
        def exec = new WrExecutor(session: session)

        File file = File.createTempFile("client", "token")
        file.deleteOnExit()
        BufferedWriter bw = new BufferedWriter(new FileWriter(file))
        String token = 'foo'
        bw.write(token)
        bw.close()

        File pem = File.createTempFile("myca", "pem")
        pem.deleteOnExit()
        bw = new BufferedWriter(new FileWriter(pem))
        String cert = '''-----BEGIN CERTIFICATE-----
MIIDDzCCAfegAwIBAgIRAJxTVSwVWoieLqbOChbcIdMwDQYJKoZIhvcNAQELBQAw
FTETMBEGA1UEChMKd3IgbWFuYWdlcjAeFw0xOTAyMTgxNTQ5NDZaFw0yMDAyMTgx
NTQ5NDZaMBUxEzARBgNVBAoTCndyIG1hbmFnZXIwggEiMA0GCSqGSIb3DQEBAQUA
A4IBDwAwggEKAoIBAQDD1zH/0QxBZhgRVEuUoP/oTvS5aOzQ7zvUhVdF7rfva4Dd
HRtXIujs2C2HmYvt3AQK6TIrscjwx3UKbD6ySQVvsS/64SRelj3Iyef48+9WHGA8
zmvVsZ+Op8LrB8gGRn+wuL51a6wXZOSM4KBachltraVRulYKCFS5h8IQIsCyZTIU
2NsDehNbT4rK6pC/kxVa9cO0+lR1NSgB4Z6h/fPjCSdOti32RTEONaZDB2LdKzd7
FaS5BfR8uWuo3A9qTxmaqgBQ3MmzI8LLvE7qOINxHqs6E6pb3G1Dzc9PhZhp4c88
A+dwP+TaFQW8yJS0Y3iWv9N/KJ9YSkgqLS48xqPPAgMBAAGjWjBYMA4GA1UdDwEB
/wQEAwICpDATBgNVHSUEDDAKBggrBgEFBQcDATAPBgNVHRMBAf8EBTADAQH/MCAG
A1UdEQQZMBeCCWxvY2FsaG9zdIcEAAAAAIcEfwAAATANBgkqhkiG9w0BAQsFAAOC
AQEAAuSMH9X1o2SJRhV+fxTWCMWu6uzNN3IZyDSkX2PjgsrQSZPLmhgbrtNE0Rfr
hjAXKPvARLUVq0cVcIYxbw/fDVOpLYvpbPpvs82d8EDWsSeHW1iP+TNk2RdU40WA
17b0zhJ1/95bTSSNQute1vnh+oqrEAEj50W9pqQd0gAiNcM1SV60VgXl3mYvIwl/
wNJ/UzQ1VKowqsbZVloElPYxs/gIDThZ4/H+xkJP2H5Igq9UZBCSUkQllHkQfhs+
MuStf8afk6GJAiH5s+N2YM/S8VWiil7oC+mq2H9PEjFKrznw95MnK5zvAi+REAYu
BNSdnfiH4cGF4cO6VXWtTfF5+Q==
-----END CERTIFICATE-----
'''
        bw.write(cert)
        bw.close()

        when:
        def result = exec.getToken()
        then:
        exec.getSession() >> session
        session.getConfigAttribute('executor.wr.tokenpath',_) >> file.path
        result == token

        when:
        def client = exec.getClient()
        then:
        exec.getSession() >> session
        session.getConfigAttribute('executor.wr.tokenpath',_) >> file.path
        session.getConfigAttribute('executor.wr.cacertpath',_) >> pem.path
        client instanceof WrRestApi
    }

}
