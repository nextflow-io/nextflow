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
 */

package nextflow.k8s.client
import javax.net.ssl.KeyManager
import java.nio.file.Files

import spock.lang.Specification
import test.TestHelper
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigDiscoveryTest extends Specification {

    def 'should read config from file' () {

        given:
        final CERT_DATA = "d29ybGQgaGVsbG8="
        final CLIENT_CERT = "aGVsbG8gd29ybGQ="
        final CLIENT_KEY = "Y2lhbyBtaWFv"

        def CONFIG = TestHelper.createInMemTempFile('config')
        CONFIG.text = """
            apiVersion: v1
            clusters:
            - cluster:
                insecure-skip-tls-verify: true
                server: https://localhost:6443
                certificate-authority-data: $CERT_DATA
              name: docker-for-desktop-cluster
            contexts:
            - context:
                cluster: docker-for-desktop-cluster
                user: docker-for-desktop
              name: docker-for-desktop
            current-context: docker-for-desktop
            kind: Config
            preferences: {}
            users:
            - name: docker-for-desktop
              user:
                client-certificate-data: $CLIENT_CERT
                client-key-data: $CLIENT_KEY
            """
            .stripIndent()

        def discovery = Spy(ConfigDiscovery)
        def KEY_MANAGERS = [] as KeyManager[]

        when:
        def config = discovery.fromKubeConfig(CONFIG, null, null, null)
        then:
        0 * discovery.discoverAuthToken(_, 'default',null) >> 'secret-token'
        1 * discovery.createKeyManagers(CLIENT_CERT.decodeBase64(), CLIENT_KEY.decodeBase64()) >>  KEY_MANAGERS
        config.server == 'https://localhost:6443'
        config.token == null
        config.namespace == 'default'
        config.serviceAccount == 'default'
        config.clientCert == CLIENT_CERT.decodeBase64()
        config.clientKey == CLIENT_KEY.decodeBase64()
        config.sslCert == CERT_DATA.decodeBase64()
        config.keyManagers.is( KEY_MANAGERS )
        !config.verifySsl
        !config.isFromCluster

    }

    def 'should read config from file with provided namespace' () {

        given:
        final CERT_DATA = "d29ybGQgaGVsbG8="
        final CLIENT_CERT = "aGVsbG8gd29ybGQ="
        final CLIENT_KEY = "Y2lhbyBtaWFv"

        def CONFIG = TestHelper.createInMemTempFile('config')
        CONFIG.text = """
            apiVersion: v1
            clusters:
            - cluster:
                insecure-skip-tls-verify: true
                server: https://localhost:6443
                certificate-authority-data: $CERT_DATA
              name: docker-for-desktop-cluster
            contexts:
            - context:
                cluster: docker-for-desktop-cluster
                user: docker-for-desktop
              name: docker-for-desktop
            current-context: docker-for-desktop
            kind: Config
            preferences: {}
            users:
            - name: docker-for-desktop
              user:
                client-certificate-data: $CLIENT_CERT
                client-key-data: $CLIENT_KEY
            """
                .stripIndent()

        def discovery = Spy(ConfigDiscovery)
        def KEY_MANAGERS = [] as KeyManager[]

        when:
        def config = discovery.fromKubeConfig(CONFIG, 'docker-for-desktop', 'ns1', 'sa2')
        then:
        0 * discovery.discoverAuthToken('docker-for-desktop','ns1','sa2') >> 'secret-token'
        1 * discovery.createKeyManagers(CLIENT_CERT.decodeBase64(), CLIENT_KEY.decodeBase64()) >>  KEY_MANAGERS
        config.server == 'https://localhost:6443'
        config.token == null
        config.namespace == 'ns1'
        config.serviceAccount == 'sa2'
        config.clientCert == CLIENT_CERT.decodeBase64()
        config.clientKey == CLIENT_KEY.decodeBase64()
        config.sslCert == CERT_DATA.decodeBase64()
        config.keyManagers.is( KEY_MANAGERS )
        !config.verifySsl
        !config.isFromCluster

    }

    def 'should read config from file with cert files' () {

        given:
        def folder = Files.createTempDirectory(null)
        def CA_FILE = folder.resolve('ca'); CA_FILE.text = 'ca-content'
        def CLIENT_CERT_FILE = folder.resolve('client-cert'); CLIENT_CERT_FILE.text = 'client-cert-content'
        def CLIENT_KEY_FILE = folder.resolve('client-key'); CLIENT_KEY_FILE.text = 'client-key-content'
        def CONFIG = folder.resolve('config')
        def KEY_MANAGERS = [] as KeyManager[]

        CONFIG.text = """
            apiVersion: v1
            clusters:
            - cluster:
                insecure-skip-tls-verify: true
                server: https://localhost:6443
                certificate-authority: $CA_FILE
              name: docker-for-desktop-cluster
            contexts:
            - context:
                cluster: docker-for-desktop-cluster
                user: docker-for-desktop
              name: docker-for-desktop
            current-context: docker-for-desktop
            kind: Config
            preferences: {}
            users:
            - name: docker-for-desktop
              user:
                client-certificate: $CLIENT_CERT_FILE
                client-key: $CLIENT_KEY_FILE
            """
                .stripIndent()

        def discovery = Spy(ConfigDiscovery)

        when:
        def config = discovery.fromKubeConfig(CONFIG, null, null, null)
        then:
        1 * discovery.createKeyManagers( CLIENT_CERT_FILE.bytes, CLIENT_KEY_FILE.bytes ) >> KEY_MANAGERS
        config.server == 'https://localhost:6443'
        config.token == null
        config.namespace == 'default'
        config.serviceAccount == 'default'
        config.clientCert == CLIENT_CERT_FILE.bytes
        config.clientKey == CLIENT_KEY_FILE.bytes
        config.sslCert == CA_FILE.bytes
        config.keyManagers.is( KEY_MANAGERS )
        !config.verifySsl
        !config.isFromCluster

        cleanup:
        folder?.deleteDir()
    }

    def 'should read config and use token' () {

        given:
        def folder = Files.createTempDirectory(null)
        def CONFIG = folder.resolve('config')

        CONFIG.text = """
            apiVersion: v1
            clusters:
            - cluster:
                insecure-skip-tls-verify: true
                server: https://localhost:6443
              name: docker-for-desktop-cluster
            contexts:
            - context:
                cluster: docker-for-desktop-cluster
                user: docker-for-desktop
              name: docker-for-desktop
            current-context: docker-for-desktop
            kind: Config
            preferences: {}
            users:
            - name: docker-for-desktop
              user:
                token: 90s090s98s7f8s
            """
                .stripIndent()

        def discovery = Spy(ConfigDiscovery)

        when:
        def config = discovery.fromKubeConfig(CONFIG, null, null, null)
        then:
        0 * discovery.discoverAuthToken(_,_,_) >> 'secret-token'
        0 * discovery.createKeyManagers( _, _ ) >> null
        config.server == 'https://localhost:6443'
        config.token == '90s090s98s7f8s'
        config.namespace == 'default'
        config.serviceAccount == 'default'
        !config.verifySsl
        !config.isFromCluster

        cleanup:
        folder?.deleteDir()
    }

    def 'should read config and discover token' () {

        given:
        def folder = Files.createTempDirectory(null)
        def CONFIG = folder.resolve('config')

        CONFIG.text = """
            apiVersion: v1
            clusters:
            - cluster:
                insecure-skip-tls-verify: true
                server: https://localhost:6443
              name: docker-for-desktop-cluster
            contexts:
            - context:
                cluster: docker-for-desktop-cluster
                user: docker-for-desktop
              name: docker-for-desktop
            current-context: docker-for-desktop
            kind: Config
            preferences: {}
            users:
            - name: docker-for-desktop
              user:
                foo: bar
            """
                .stripIndent()

        def discovery = Spy(ConfigDiscovery)

        when:
        def config = discovery.fromKubeConfig(CONFIG, null, null, null)
        then:
        1 * discovery.discoverAuthToken(_, _, _) >> 'secret-token'
        0 * discovery.createKeyManagers( _, _ ) >> null
        config.server == 'https://localhost:6443'
        config.token == 'secret-token'
        config.namespace == 'default'
        config.serviceAccount == 'default'
        !config.verifySsl
        !config.isFromCluster

        cleanup:
        folder?.deleteDir()
    }

    def 'should read config from given context' () {

        given:
        def folder = Files.createTempDirectory(null)
        folder.resolve('fake-cert-file').text = 'fake-cert-content'
        folder.resolve('fake-key-file').text = 'fake-key-content'
        folder.resolve('fake-ca-file').text = 'fake-ca-content'

        def CONFIG = folder.resolve('config')
        CONFIG.text = '''
                apiVersion: v1
                clusters:
                - cluster:
                    certificate-authority: fake-ca-file
                    server: https://1.2.3.4
                  name: development
                - cluster:
                    insecure-skip-tls-verify: true
                    server: https://5.6.7.8
                  name: scratch
                contexts:
                - context:
                    cluster: development
                    namespace: frontend
                    user: developer
                  name: dev-frontend
                - context:
                    cluster: development
                    namespace: storage
                    user: developer
                  name: dev-storage
                - context:
                    cluster: scratch
                    namespace: default
                    user: experimenter
                  name: exp-scratch
                current-context: ""
                kind: Config
                preferences: {}
                users:
                - name: developer
                  user:
                    client-certificate: fake-cert-file
                - name: experimenter
                  user:
                    password: some-password
                    username: exp        
        '''.stripIndent()

        when:
        def cfg1 = new ConfigDiscovery().fromKubeConfig(CONFIG, 'dev-frontend', null, null)
        then:
        cfg1.server == 'https://1.2.3.4'
        cfg1.sslCert == 'fake-ca-content'.bytes
        cfg1.isVerifySsl()
        cfg1.namespace == 'frontend'
        cfg1.serviceAccount == 'default'
        cfg1.clientCert == 'fake-cert-content'.bytes

        when:
        def cfg2 = new ConfigDiscovery().fromKubeConfig(CONFIG, 'dev-storage', null, null)
        then:
        cfg2.server == 'https://1.2.3.4'
        cfg2.sslCert == 'fake-ca-content'.bytes
        cfg2.isVerifySsl()
        cfg2.namespace == 'storage'
        cfg2.serviceAccount == 'default'
        cfg2.clientCert == 'fake-cert-content'.bytes

        when:
        def cfg3 = new ConfigDiscovery().fromKubeConfig(CONFIG, 'exp-scratch', null, null)
        then:
        cfg3.server == 'https://5.6.7.8'
        cfg3.sslCert == null
        !cfg3.isVerifySsl()
        cfg3.namespace == 'default'
        cfg3.serviceAccount == 'default'

        when:
        new ConfigDiscovery().fromKubeConfig(CONFIG, 'foo', null, null)
        then:
        thrown(IllegalArgumentException)

        true
        cleanup:
        folder.deleteDir()
    }

    def 'should load from cluster env' () {
        given:
        def CERT_FILE = TestHelper.createInMemTempFile('ca'); CERT_FILE.text = 'ca-content'
        def TOKEN_FILE = TestHelper.createInMemTempFile('token'); TOKEN_FILE.text = 'my-token'
        def NAMESPACE_FILE = TestHelper.createInMemTempFile('namespace'); NAMESPACE_FILE.text = 'foo-namespace'

        def discovery = Spy(ConfigDiscovery)

        when:
        def env = [ KUBERNETES_SERVICE_HOST: 'foo.com', KUBERNETES_SERVICE_PORT: '4343' ]
        def config = discovery.fromCluster(env, null, null)
        then:
        1 * discovery.path('/var/run/secrets/kubernetes.io/serviceaccount/ca.crt') >> CERT_FILE
        1 * discovery.path('/var/run/secrets/kubernetes.io/serviceaccount/token') >> TOKEN_FILE
        1 * discovery.path('/var/run/secrets/kubernetes.io/serviceaccount/namespace') >> NAMESPACE_FILE
        0 * discovery.createKeyManagers(_,_) >> null
        and:
        config.server == 'foo.com:4343'
        config.namespace == 'foo-namespace'
        config.token == 'my-token'
        config.sslCert == CERT_FILE.text.bytes
        config.isFromCluster

        when:
        env = [ KUBERNETES_SERVICE_HOST: 'https://host.com' ]
        config = discovery.fromCluster(env, 'my-namespace', null)
        then:
        1 * discovery.path('/var/run/secrets/kubernetes.io/serviceaccount/ca.crt') >> CERT_FILE
        1 * discovery.path('/var/run/secrets/kubernetes.io/serviceaccount/token') >> TOKEN_FILE
        1 * discovery.path('/var/run/secrets/kubernetes.io/serviceaccount/namespace') >> NAMESPACE_FILE
        and:
        config.server == 'https://host.com'
        config.namespace == 'my-namespace'
    }

    def 'should create  key managers' () {
        given:
        final CERT = 'LS0tLS1CRUdJTiBDRVJUSUZJQ0FURS0tLS0tCk1JSUNTRENDQVRDZ0F3SUJBZ0lJRlFsM1l2Y2k1TWN3RFFZSktvWklodmNOQVFFTEJRQXdGVEVUTUJFR0ExVUUKQXhNS2EzVmlaWEp1WlhSbGN6QWVGdzB4T0RBeE1UTXlNREEyTlRaYUZ3MHhPVEF4TVRNeU1EQTJOVFphTURNeApGREFTQmdOVkJBb1RDMFJ2WTJ0bGNpQkpibU11TVJzd0dRWURWUVFERXhKa2IyTnJaWEl0Wm05eUxXUmxjMnQwCmIzQXdnWjh3RFFZSktvWklodmNOQVFFQkJRQURnWTBBTUlHSkFvR0JBUEtYT0ZsV2t2THIzb29ETGNFOElyME0KTzNBMHZqQlVvUzZ0bUdBbFRYYTd0QWQwM3BTMXNJNit0WVRwVlU2YXR6ZU9vU0VrOWhmaWxBdVNYdG1hSHZCUAp1czFEcG1LZEZRMWI3OFRkSnQ4OGV3c3BRajFxYUwvQldHeitMUzUrRHUrNUJuUGtmZlhDS1UxQTdUc2tZamJyClhxeDhlN2FWZURWTmFjZXc0Z0RqQWdNQkFBR2pBakFBTUEwR0NTcUdTSWIzRFFFQkN3VUFBNElCQVFBMXVtVlAKR29EZTVCRXJrb21qWXdITXhiTTd4UStibTYrUDE1T0pINUo0UGNQeU11d25ocC9ORVp1NnpsTTZSUUo3SUNKQgpHWTRBMnFKVmJsWUkwQkJzRkF1TXMreTAyazdVVVVoK0NRYVd0SXhBcFNmbkQ4dUVXQ0g5VE1ZNGdLbTZjTDhVCk1OVVl1RnpUQ2hmTS96RjdUMXVaZWxJYXNrYXFaWSt3a3hxa3YyRUQxQ2F5MDUxSXRWRXZVbDIvSVZyVHdrT20KZ25nL3Q4L2RkeDhpOUkzTFJrMTlTaERKdXlQZ1NrTTZRSWlSd09mRHk4V0ZFaURpd0hBS0ErSEZhTGhOOFJTMwpieDUvdEhEN01id0FpdnorNTU4YUFEQjNEd1ZpekthM2d5Wm4yUzRjUGFqZnNwODFqRkNIQS9QekdQdTU2MzJwCkxRN0gyRW1aYmJuUHFYTFgKLS0tLS1FTkQgQ0VSVElGSUNBVEUtLS0tLQo='.decodeBase64()
        final KEY = 'LS0tLS1CRUdJTiBSU0EgUFJJVkFURSBLRVktLS0tLQpNSUlDWGdJQkFBS0JnUUR5bHpoWlZwTHk2OTZLQXkzQlBDSzlERHR3Tkw0d1ZLRXVyWmhnSlUxMnU3UUhkTjZVCnRiQ092cldFNlZWT21yYzNqcUVoSlBZWDRwUUxrbDdabWg3d1Q3ck5RNlppblJVTlcrL0UzU2JmUEhzTEtVSTkKYW1pL3dWaHMvaTB1Zmc3dnVRWno1SDMxd2lsTlFPMDdKR0kyNjE2c2ZIdTJsWGcxVFduSHNPSUE0d0lEQVFBQgpBb0dBYWRUOCtVU2lvU1d6bFVRanZ1eHNQMHRKMXY2N2hqdzFnVGFzaGkxZjZRK2tUNmgxdml5eGxPU3dMZ2JaCmQ0eFpwL3dxWVZwTm5rZnp6RVNUNnB5cEo5WTEwdHY1cFpSWG9HbG1NT2tIZSswUW45N0c5ZDRzL2JCV3lmYXYKRzhRTC9tZFN6Vy85YUdrSkpiNWU0VDlsSURvRDNFVDgwYUFWbzl2V0NPVUxsdWtDUVFEK0hINU5ucVBuSTdnTApWOUJKZzlRRVBwUTVYa2traW8rejZ2YkRHQU5rR1VPV1dmRURKUHE2Q2JBb1dqeWh1Qy9KS1dYRWs4Rkt0M1Y2CkhVNllYeVpGQWtFQTlHVE9XOFM4KzVNNHE5R3lNeURxN1ZkVHA3M2daeSsvNjVQam5hNlpDUnhTZklxL2xKUVoKY2F6MkhGYVRzRFdLbkdhWGNxTmdBVXNEODNyWTlzM3hCd0pCQUt5Vjc1YUtPMm0rRWI3cWVsV2p5bmpEZytwZQp4akNpUnkxOFZQSjJPYjlmaFU3MWNVS2dlQVdvbE5NalRuREw1dkNxUkNzNTZ4cnk5VC9sN2I2QlNUMENRUURnCjRoV2xDZTdnQzhOZEQzTkxhdUhpRGJZenB4dmp0Mk9Ca2E4ai9ISmptTVVxUnI0dEtPNFUxUlFPVlhoRzc2MmgKWnlHNjRpeklZOCs1N3ZQUWZ3Wm5Ba0VBdW9RWW1lUi90UWhIakhRNFlhZGRHbkNBQ2hZZ29ObEFzSGhGTElxVQo1ZTZaMXN2Q3VKU285TDVVRCtrclFUYWlGU01pRHZwZlJyVE1ZKzZ5Q0tTajd3PT0KLS0tLS1FTkQgUlNBIFBSSVZBVEUgS0VZLS0tLS0K'.decodeBase64()
        final discovery = new ConfigDiscovery()
        when:
        def managers = discovery.createKeyManagers(CERT, KEY)
        then:
        managers.size()==1
    }

    def 'should create key managers from an EC client key' () {
        given:
        final CERT = 'LS0tLS1CRUdJTiBDRVJUSUZJQ0FURS0tLS0tCk1JSUNTRENDQVRDZ0F3SUJBZ0lJRlFsM1l2Y2k1TWN3RFFZSktvWklodmNOQVFFTEJRQXdGVEVUTUJFR0ExVUUKQXhNS2EzVmlaWEp1WlhSbGN6QWVGdzB4T0RBeE1UTXlNREEyTlRaYUZ3MHhPVEF4TVRNeU1EQTJOVFphTURNeApGREFTQmdOVkJBb1RDMFJ2WTJ0bGNpQkpibU11TVJzd0dRWURWUVFERXhKa2IyTnJaWEl0Wm05eUxXUmxjMnQwCmIzQXdnWjh3RFFZSktvWklodmNOQVFFQkJRQURnWTBBTUlHSkFvR0JBUEtYT0ZsV2t2THIzb29ETGNFOElyME0KTzNBMHZqQlVvUzZ0bUdBbFRYYTd0QWQwM3BTMXNJNit0WVRwVlU2YXR6ZU9vU0VrOWhmaWxBdVNYdG1hSHZCUAp1czFEcG1LZEZRMWI3OFRkSnQ4OGV3c3BRajFxYUwvQldHeitMUzUrRHUrNUJuUGtmZlhDS1UxQTdUc2tZamJyClhxeDhlN2FWZURWTmFjZXc0Z0RqQWdNQkFBR2pBakFBTUEwR0NTcUdTSWIzRFFFQkN3VUFBNElCQVFBMXVtVlAKR29EZTVCRXJrb21qWXdITXhiTTd4UStibTYrUDE1T0pINUo0UGNQeU11d25ocC9ORVp1NnpsTTZSUUo3SUNKQgpHWTRBMnFKVmJsWUkwQkJzRkF1TXMreTAyazdVVVVoK0NRYVd0SXhBcFNmbkQ4dUVXQ0g5VE1ZNGdLbTZjTDhVCk1OVVl1RnpUQ2hmTS96RjdUMXVaZWxJYXNrYXFaWSt3a3hxa3YyRUQxQ2F5MDUxSXRWRXZVbDIvSVZyVHdrT20KZ25nL3Q4L2RkeDhpOUkzTFJrMTlTaERKdXlQZ1NrTTZRSWlSd09mRHk4V0ZFaURpd0hBS0ErSEZhTGhOOFJTMwpieDUvdEhEN01id0FpdnorNTU4YUFEQjNEd1ZpekthM2d5Wm4yUzRjUGFqZnNwODFqRkNIQS9QekdQdTU2MzJwCkxRN0gyRW1aYmJuUHFYTFgKLS0tLS1FTkQgQ0VSVElGSUNBVEUtLS0tLQo='.decodeBase64()
        final KEY = 'LS0tLS1CRUdJTiBQUklWQVRFIEtFWS0tLS0tCk1JR0hBZ0VBTUJNR0J5cUdTTTQ5QWdFR0NDcUdTTTQ5QXdFSEJHMHdhd0lCQVFRZ21aZFZ3NmJRU0w1T1l5RjQKbzJ4V0hUQ05BSW1hRTkycGd2dGMzK2Z2UDVxaFJBTkNBQVJSd0RpUVptTUNqcWxvbFBzRTdiZjgwWjhrZkRXTworS2U4NUdVSll2MlBubWVxbDhkYjdwcmFlMHFPQUJaaXR2Mmh2SmJFeFdsUFR0MS9CYTNMK1B5NAotLS0tLUVORCBQUklWQVRFIEtFWS0tLS0tCg=='.decodeBase64()
        final discovery = new ConfigDiscovery()
        when:
        def managers = discovery.createKeyManagers(CERT, KEY)
        then:
        managers.size()==1
    }

    def 'should create key managers from an EC-encrypted client key' () {
        given:
        final CERT = 'LS0tLS1CRUdJTiBDRVJUSUZJQ0FURS0tLS0tCk1JSUJrVENDQVRlZ0F3SUJBZ0lJSGw1Zmx0UmRTdDB3Q2dZSUtvWkl6ajBFQXdJd0l6RWhNQjhHQTFVRUF3d1kKYXpOekxXTnNhV1Z1ZEMxallVQXhOekUwTnpRM09UVTNNQjRYRFRJME1EVXdNekUwTlRJek4xb1hEVEkxTURVdwpNekUwTlRJek4xb3dNREVYTUJVR0ExVUVDaE1PYzNsemRHVnRPbTFoYzNSbGNuTXhGVEFUQmdOVkJBTVRESE41CmMzUmxiVHBoWkcxcGJqQlpNQk1HQnlxR1NNNDlBZ0VHQ0NxR1NNNDlBd0VIQTBJQUJQN1Q5RHVvUlllLzBlUkwKUmNHV2RoYnl2Q3BucXlsSVIyaUwxdGkwc1hVdEpZZjUrVXhIOWFBMjdzY2FSYW1qbjdnTTFrKzZNaVk5cm15OApyRmdoWm1xalNEQkdNQTRHQTFVZER3RUIvd1FFQXdJRm9EQVRCZ05WSFNVRUREQUtCZ2dyQmdFRkJRY0RBakFmCkJnTlZIU01FR0RBV2dCU0NhdXFoQVEvWEdoaFRtaFBoY21vRVdOeWluakFLQmdncWhrak9QUVFEQWdOSUFEQkYKQWlCZzRaNmlWeFV3Mk5uMHBQTG02VlovUGttQnVuTDEwZG50dEg3UVdIcklCd0loQU0vTDhVMGxQN0IyeFEyZwpsZjlhNHNhbzJ1bE5ONnQvQ0dibzlxTlo1QzZHCi0tLS0tRU5EIENFUlRJRklDQVRFLS0tLS0KLS0tLS1CRUdJTiBDRVJUSUZJQ0FURS0tLS0tCk1JSUJkekNDQVIyZ0F3SUJBZ0lCQURBS0JnZ3Foa2pPUFFRREFqQWpNU0V3SHdZRFZRUUREQmhyTTNNdFkyeHAKWlc1MExXTmhRREUzTVRRM05EYzVOVGN3SGhjTk1qUXdOVEF6TVRRMU1qTTNXaGNOTXpRd05UQXhNVFExTWpNMwpXakFqTVNFd0h3WURWUVFEREJock0zTXRZMnhwWlc1MExXTmhRREUzTVRRM05EYzVOVGN3V1RBVEJnY3Foa2pPClBRSUJCZ2dxaGtqT1BRTUJCd05DQUFRZHFYVHdIQS9mVjRKZGdYa2FubXB1OVE0QStwUGRGaXZGdytiUmVhdEYKUXVOUTBKWndIbzlaa2ltb2lEUU5qb2h0TWdHckdtTVlsTTZuaXM4ZVFvM3RvMEl3UURBT0JnTlZIUThCQWY4RQpCQU1DQXFRd0R3WURWUjBUQVFIL0JBVXdBd0VCL3pBZEJnTlZIUTRFRmdRVWdtcnFvUUVQMXhvWVU1b1Q0WEpxCkJGamNvcDR3Q2dZSUtvWkl6ajBFQXdJRFNBQXdSUUloQUlvb2ZmNzdvb1VYS2hmNVo3aVRzdExhOTVwU2VaRmUKRHZjMXdFQXVEa3NTQWlBNzJQajJxNnpBclhpYkpUa0s2RTBHTEtVODdhTHhHc3BmS29uVVJnalI2Zz09Ci0tLS0tRU5EIENFUlRJRklDQVRFLS0tLS0K'.decodeBase64()
        final KEY = 'LS0tLS1CRUdJTiBFQyBQUklWQVRFIEtFWS0tLS0tCk1IY0NBUUVFSUNvQTNvRHkzN3NXdmszM3JGRGtRdlZ1Wkh1cCt1Uk40V3RqbUlPR1c4cHBvQW9HQ0NxR1NNNDkKQXdFSG9VUURRZ0FFL3RQME82aEZoNy9SNUV0RndaWjJGdks4S21lcktVaEhhSXZXMkxTeGRTMGxoL241VEVmMQpvRGJ1eHhwRnFhT2Z1QXpXVDdveUpqMnViTHlzV0NGbWFnPT0KLS0tLS1FTkQgRUMgUFJJVkFURSBLRVktLS0tLQo='.decodeBase64()
        final discovery = new ConfigDiscovery()
        when:
        def managers = discovery.createKeyManagers(CERT, KEY)
        then:
        managers.size()==1
    }
}
