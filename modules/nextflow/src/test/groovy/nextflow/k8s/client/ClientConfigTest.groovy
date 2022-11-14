/*
 * Copyright 2020-2022, Seqera Labs
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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

import java.nio.file.Files

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ClientConfigTest extends Specification {

    def 'should stringify a config' () {

        when:
        final CERT = 'LS0tLS1CRUdJTiBDRVJUSUZJQ0FURS0tLS0tCk1JSUNTRENDQVRDZ0F3SUJBZ0lJRlFsM1l2Y2k1TWN3RFFZSktvWklodmNOQVFFTEJRQXdGVEVUTUJFR0ExVUUKQXhNS2EzVmlaWEp1WlhSbGN6QWVGdzB4T0RBeE1UTXlNREEyTlRaYUZ3MHhPVEF4TVRNeU1EQTJOVFphTURNeApGREFTQmdOVkJBb1RDMFJ2WTJ0bGNpQkpibU11TVJzd0dRWURWUVFERXhKa2IyTnJaWEl0Wm05eUxXUmxjMnQwCmIzQXdnWjh3RFFZSktvWklodmNOQVFFQkJRQURnWTBBTUlHSkFvR0JBUEtYT0ZsV2t2THIzb29ETGNFOElyME0KTzNBMHZqQlVvUzZ0bUdBbFRYYTd0QWQwM3BTMXNJNit0WVRwVlU2YXR6ZU9vU0VrOWhmaWxBdVNYdG1hSHZCUAp1czFEcG1LZEZRMWI3OFRkSnQ4OGV3c3BRajFxYUwvQldHeitMUzUrRHUrNUJuUGtmZlhDS1UxQTdUc2tZamJyClhxeDhlN2FWZURWTmFjZXc0Z0RqQWdNQkFBR2pBakFBTUEwR0NTcUdTSWIzRFFFQkN3VUFBNElCQVFBMXVtVlAKR29EZTVCRXJrb21qWXdITXhiTTd4UStibTYrUDE1T0pINUo0UGNQeU11d25ocC9ORVp1NnpsTTZSUUo3SUNKQgpHWTRBMnFKVmJsWUkwQkJzRkF1TXMreTAyazdVVVVoK0NRYVd0SXhBcFNmbkQ4dUVXQ0g5VE1ZNGdLbTZjTDhVCk1OVVl1RnpUQ2hmTS96RjdUMXVaZWxJYXNrYXFaWSt3a3hxa3YyRUQxQ2F5MDUxSXRWRXZVbDIvSVZyVHdrT20KZ25nL3Q4L2RkeDhpOUkzTFJrMTlTaERKdXlQZ1NrTTZRSWlSd09mRHk4V0ZFaURpd0hBS0ErSEZhTGhOOFJTMwpieDUvdEhEN01id0FpdnorNTU4YUFEQjNEd1ZpekthM2d5Wm4yUzRjUGFqZnNwODFqRkNIQS9QekdQdTU2MzJwCkxRN0gyRW1aYmJuUHFYTFgKLS0tLS1FTkQgQ0VSVElGSUNBVEUtLS0tLQo='.decodeBase64()
        final KEY = 'LS0tLS1CRUdJTiBSU0EgUFJJVkFURSBLRVktLS0tLQpNSUlDWGdJQkFBS0JnUUR5bHpoWlZwTHk2OTZLQXkzQlBDSzlERHR3Tkw0d1ZLRXVyWmhnSlUxMnU3UUhkTjZVCnRiQ092cldFNlZWT21yYzNqcUVoSlBZWDRwUUxrbDdabWg3d1Q3ck5RNlppblJVTlcrL0UzU2JmUEhzTEtVSTkKYW1pL3dWaHMvaTB1Zmc3dnVRWno1SDMxd2lsTlFPMDdKR0kyNjE2c2ZIdTJsWGcxVFduSHNPSUE0d0lEQVFBQgpBb0dBYWRUOCtVU2lvU1d6bFVRanZ1eHNQMHRKMXY2N2hqdzFnVGFzaGkxZjZRK2tUNmgxdml5eGxPU3dMZ2JaCmQ0eFpwL3dxWVZwTm5rZnp6RVNUNnB5cEo5WTEwdHY1cFpSWG9HbG1NT2tIZSswUW45N0c5ZDRzL2JCV3lmYXYKRzhRTC9tZFN6Vy85YUdrSkpiNWU0VDlsSURvRDNFVDgwYUFWbzl2V0NPVUxsdWtDUVFEK0hINU5ucVBuSTdnTApWOUJKZzlRRVBwUTVYa2traW8rejZ2YkRHQU5rR1VPV1dmRURKUHE2Q2JBb1dqeWh1Qy9KS1dYRWs4Rkt0M1Y2CkhVNllYeVpGQWtFQTlHVE9XOFM4KzVNNHE5R3lNeURxN1ZkVHA3M2daeSsvNjVQam5hNlpDUnhTZklxL2xKUVoKY2F6MkhGYVRzRFdLbkdhWGNxTmdBVXNEODNyWTlzM3hCd0pCQUt5Vjc1YUtPMm0rRWI3cWVsV2p5bmpEZytwZQp4akNpUnkxOFZQSjJPYjlmaFU3MWNVS2dlQVdvbE5NalRuREw1dkNxUkNzNTZ4cnk5VC9sN2I2QlNUMENRUURnCjRoV2xDZTdnQzhOZEQzTkxhdUhpRGJZenB4dmp0Mk9Ca2E4ai9ISmptTVVxUnI0dEtPNFUxUlFPVlhoRzc2MmgKWnlHNjRpeklZOCs1N3ZQUWZ3Wm5Ba0VBdW9RWW1lUi90UWhIakhRNFlhZGRHbkNBQ2hZZ29ObEFzSGhGTElxVQo1ZTZaMXN2Q3VKU285TDVVRCtrclFUYWlGU01pRHZwZlJyVE1ZKzZ5Q0tTajd3PT0KLS0tLS1FTkQgUlNBIFBSSVZBVEUgS0VZLS0tLS0K'.decodeBase64()

        def config = new ClientConfig()
        config.clientCert = CERT
        config.clientKey = KEY
        config.sslCert = CERT
        println config.toString()

        then:
        noExceptionThrown()

    }

    def 'should create a client config from a map' () {

        given:
        def MAP = [
                server:'foo.com',
                token: 'blah-blah',
                namespace: 'my-namespace',
                verifySsl: true,
                sslCert: 'fizzbuzz'.bytes.encodeBase64().toString(),
                clientCert: 'hello'.bytes.encodeBase64().toString(),
                clientKey: 'world'.bytes.encodeBase64().toString() ]

        when:
        def result = ClientConfig.fromNextflowConfig(MAP, null, null)

        then:
        result.server == 'foo.com'
        result.token == 'blah-blah'
        result.namespace == 'my-namespace'
        result.serviceAccount == 'default'
        result.verifySsl
        result.clientCert == 'hello'.bytes
        result.clientKey == 'world'.bytes
        result.sslCert == 'fizzbuzz'.bytes

        when:
        result = ClientConfig.fromNextflowConfig(MAP, 'ns1', 'sa2')
        then:
        result.server == 'foo.com'
        result.token == 'blah-blah'
        result.namespace == 'ns1'
        result.serviceAccount == 'sa2'
        result.verifySsl
        result.clientCert == 'hello'.bytes
        result.clientKey == 'world'.bytes
        result.sslCert == 'fizzbuzz'.bytes
    }

    def 'should create a client config from a map with files' () {

        given:
        def folder = Files.createTempDirectory('test')
        def file1 = folder.resolve('file1')
        def file2 = folder.resolve('file2')
        def file3 = folder.resolve('file3')
        file1.text = 'fizzbuzz'.bytes.encodeBase64().toString()
        file2.text = 'hello'.bytes.encodeBase64().toString()
        file3.text = 'world'.bytes.encodeBase64().toString()

        def MAP = [
                server:'foo.com',
                token: 'blah-blah',
                namespace: 'my-namespace',
                verifySsl: false,
                sslCertFile: file1,
                clientCertFile: file2,
                clientKeyFile: file3 ]

        when:
        def result = ClientConfig.fromNextflowConfig(MAP, null, null)

        then:
        result.server == 'foo.com'
        result.token == 'blah-blah'
        result.namespace == 'my-namespace'
        result.serviceAccount == 'default'
        !result.verifySsl
        result.sslCert == file1.text.bytes
        result.clientCert == file2.text.bytes
        result.clientKey == file3.text.bytes

        cleanup:
        folder?.deleteDir()
    }

}
