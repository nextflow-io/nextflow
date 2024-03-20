/*
 * Copyright 2013-2023, Seqera Labs
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
import spock.lang.Unroll
import spock.util.environment.RestoreSystemProperties
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProxyHelperTest extends Specification {

    @RestoreSystemProperties
    def 'should setup proxy properties'() {

        when:
        ProxyHelper.setProxy('HTTP', [HTTP_PROXY: 'alpha.com:333'])
        then:
        System.getProperty('http.proxyHost') == 'alpha.com'
        System.getProperty('http.proxyPort') == '333'

        when:
        ProxyHelper.setProxy('http', [http_proxy: 'gamma.com:444'])
        then:
        System.getProperty('http.proxyHost') == 'gamma.com'
        System.getProperty('http.proxyPort') == '444'

        when:
        ProxyHelper.setProxy('HTTPS', [HTTPS_PROXY: 'beta.com:5466'])
        then:
        System.getProperty('https.proxyHost') == 'beta.com'
        System.getProperty('https.proxyPort') == '5466'

        when:
        ProxyHelper.setProxy('https', [https_proxy: 'zeta.com:6646'])
        then:
        System.getProperty('https.proxyHost') == 'zeta.com'
        System.getProperty('https.proxyPort') == '6646'

        when:
        ProxyHelper.setProxy('FTP', [FTP_PROXY: 'delta.com:7566'])
        then:
        System.getProperty('ftp.proxyHost') == 'delta.com'
        System.getProperty('ftp.proxyPort') == '7566'

        when:
        ProxyHelper.setProxy('ftp', [ftp_proxy: 'epsilon.com:6658'])
        then:
        System.getProperty('ftp.proxyHost') == 'epsilon.com'
        System.getProperty('ftp.proxyPort') == '6658'
    }

    @RestoreSystemProperties
    def 'should setup proxy properties and configure the network authenticator'() {

        when:
        ProxyHelper.setProxy('HTTP', [HTTP_PROXY: 'http://alphauser:alphapass@alpha.com:333'])
        PasswordAuthentication auth = Authenticator.requestPasswordAuthentication(
            'alpha.com', null, 333, 'http', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('http.proxyHost') == 'alpha.com'
        System.getProperty('http.proxyPort') == '333'
        and:
        auth.getUserName() == 'alphauser'
        auth.getPassword() == 'alphapass'.toCharArray()

        when:
        ProxyHelper.setProxy('http', [http_proxy: 'http://gammauser:gammapass@gamma.com:444'])
        auth = Authenticator.requestPasswordAuthentication(
            'gamma.com', null, 444, 'http', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('http.proxyHost') == 'gamma.com'
        System.getProperty('http.proxyPort') == '444'
        and:
        auth.getUserName() == 'gammauser'
        auth.getPassword() == 'gammapass'.toCharArray()

        when:
        ProxyHelper.setProxy('HTTPS', [HTTPS_PROXY: 'https://betauser:betapass@beta.com:5466'])
        auth = Authenticator.requestPasswordAuthentication(
            'beta.com', null, 5466, 'HTTPS', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('https.proxyHost') == 'beta.com'
        System.getProperty('https.proxyPort') == '5466'
        and:
        auth.getUserName() == 'betauser'
        auth.getPassword() == 'betapass'.toCharArray()

        when:
        ProxyHelper.setProxy('https', [https_proxy: 'https://zetauser:zetapass@zeta.com:6646'])
        auth = Authenticator.requestPasswordAuthentication(
            'zeta.com', null, 6646, 'https', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('https.proxyHost') == 'zeta.com'
        System.getProperty('https.proxyPort') == '6646'
        and:
        auth.getUserName() == 'zetauser'
        auth.getPassword() == 'zetapass'.toCharArray()

        when:
        ProxyHelper.setProxy('FTP', [FTP_PROXY: 'ftp://deltauser:deltapass@delta.com:7566'])
        auth = Authenticator.requestPasswordAuthentication(
            'delta.com', null, 7566, 'ftp', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('ftp.proxyHost') == 'delta.com'
        System.getProperty('ftp.proxyPort') == '7566'
        and:
        auth.getUserName() == 'deltauser'
        auth.getPassword() == 'deltapass'.toCharArray()

        when:
        ProxyHelper.setProxy('ftp', [ftp_proxy: 'ftp://epsilonuser:epsilonpass@epsilon.com:6658'])
        auth = Authenticator.requestPasswordAuthentication(
            'epsilon.com', null, 6658, 'ftp', null, null, null, Authenticator.RequestorType.PROXY
        )
        then:
        System.getProperty('ftp.proxyHost') == 'epsilon.com'
        System.getProperty('ftp.proxyPort') == '6658'
        and:
        auth.getUserName() == 'epsilonuser'
        auth.getPassword() == 'epsilonpass'.toCharArray()
    }

    @RestoreSystemProperties
    def 'should set no proxy property' () {

        given:
        System.properties.remove('http.nonProxyHosts')
        
        when:
        ProxyHelper.setNoProxy(ENV)
        then:
        System.getProperty('http.nonProxyHosts') == EXPECTED

        where:
        ENV                         | EXPECTED
        [:]                         | null
        [no_proxy: 'localhost' ]    | 'localhost'
        [NO_PROXY: '127.0.0.1' ]    | '127.0.0.1'
        [NO_PROXY:'localhost,127.0.0.1,.localdomain.com']  | 'localhost|127.0.0.1|.localdomain.com'

    }

    @RestoreSystemProperties
    @Unroll
    def 'should set http client timeout' () {
        when:
        ProxyHelper.setHttpClientProperties(ENV)
        then:
        System.getProperty('jdk.httpclient.keepalive.timeout') == TIMEOUT
        and:
        System.getProperty('jdk.httpclient.connectionPoolSize') == POOLSIZE

        where:
        ENV                                             | TIMEOUT   | POOLSIZE
        [:]                                             | '10'      | null
        and:
        [NXF_JDK_HTTPCLIENT_KEEPALIVE_TIMEOUT: '1']     | '1'       | null
        [NXF_JDK_HTTPCLIENT_KEEPALIVE_TIMEOUT: '100']   | '100'     | null
        and:
        [NXF_JDK_HTTPCLIENT_CONNECTIONPOOLSIZE: '0']    | '10'      | '0'
        [NXF_JDK_HTTPCLIENT_CONNECTIONPOOLSIZE: '99']   | '10'      | '99'
    }

}
