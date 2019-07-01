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

package nextflow.wr.client

import java.io.File
import java.nio.file.Paths
import java.nio.file.Path
import javax.net.ssl.HttpsURLConnection

import nextflow.processor.TaskRun
import nextflow.processor.TaskBean
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskConfig
import nextflow.wr.executor.WrFileCopyStrategy

import spock.lang.Specification
/**
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * based on K8sClientTest by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class WrRestApiTest extends Specification {

    private final certData = '''
        -----BEGIN CERTIFICATE-----
        MIIDDjCCAfagAwIBAgIQfzLCXb/maNHYDLn5PsABGzANBgkqhkiG9w0BAQsFADAV
        MRMwEQYDVQQKEwp3ciBtYW5hZ2VyMB4XDTE5MDQwNDEzNDIwNloXDTIwMDQwMzEz
        NDIwNlowFTETMBEGA1UEChMKd3IgbWFuYWdlcjCCASIwDQYJKoZIhvcNAQEBBQAD
        ggEPADCCAQoCggEBAMONmay5THFDnygpnLwwuJyno7BR1IL5myMAvBkQW5E4nLOV
        EzP7UpNUQJEuDX2mBJN1vz0kh+f+Vc4PG/7B2jUJSI9A8FlUCoSrRyLWPWtfGCom
        o41yr8YzhmWqhc+ySwDHNyXrUY/rQyjfo/tcU242GTk83QyodGE8S+Vah37tudXu
        5Q8q9HFjAd8KGak9ZPWctHueoiAp8u9fPR6m5f3aaWemDG8Fmwh9rfOJntzc/VAt
        QrtiUN7Iv7sUm5ecV2A6r8s+RKlNUxHV99q+ILsYHS7Pk1QmCAit/iJOFkosEIPX
        QIUPrQ0ALnm/9w3SrMFitwaa/fSVFPuzfl+Vjk0CAwEAAaNaMFgwDgYDVR0PAQH/
        BAQDAgKkMBMGA1UdJQQMMAoGCCsGAQUFBwMBMA8GA1UdEwEB/wQFMAMBAf8wIAYD
        VR0RBBkwF4IJbG9jYWxob3N0hwQAAAAAhwR/AAABMA0GCSqGSIb3DQEBCwUAA4IB
        AQCRgVQ+mPloPOsjZghZ/4IgHf69cFqPfjU9xudOQqNQ4YmOU/TodLUweqS2/dW6
        w10OXiOb5phaD+ss9+5MO7ZE0jpgVgPiDc8dXq6DEiDwXTgJbJtVxPUy7UxsisRC
        uBPVmoCSKmWg2J0g0+DAYZzhGgVEzEG1KQ3/ZSsiFpEsnzAl2s7DdAsSUU336uN3
        slCUCODtiYwxq2j1NT19GUHJheRAhjqd4RxsTlrQWILdR0Ohu1xmsFE0bCncn89J
        pBjTyfNvF3RL2LRSA8FD/xN/DlzQoDTlg7/mNRJ9GYfmTjnZqBlex3vnKweXudZL
        RAkgePpl2m7K9RfOJrlMCAgz
        -----END CERTIFICATE-----
        '''.stripIndent().rightTrim()

    def 'should create a client' () {
        given:
        // create a cacert file
        File capem = File.createTempFile("test", "ca.pem")
        capem.write(certData)

        String endpoint = "https://localhost:12345"
        String token = "token"

        when:
        def client = new WrRestApi(endpoint, token, capem.absolutePath)

        then:
        client.endpoint == endpoint
        client.token == token

        cleanup:
        capem.delete()
    }

    private WrRestApi createClient() {
        def capem = File.createTempFile("test","ca.pem").with {
            deleteOnExit()
            write certData
            return absolutePath
        }
        return new WrRestApi("https://localhost:12345", "token", capem)
    }

    def 'should add tasks' () {
        given:
        def client = Spy(createClient())
        def processor = new TaskProcessor()
        processor.name = 'myprocess'

        def workDir1 = Paths.get("/my/work/dir")
        Map<String,Object> entries1 = [:]
        def config1 = new TaskConfig(entries1)
        def task1 = [name: 'Hello 1', processor: processor, config: config1, workDir: workDir1] as TaskRun
        def bean1 = new TaskBean(workDir: workDir1)
        def strategy1 = new WrFileCopyStrategy(bean1, null)

        def workDir2 = 's3://bucket/work' as Path
        Map<String,Object> entries2 = [:]
        entries2 << [cpus:2]
        entries2 << [memory:"500MB"]
        entries2 << [time:"2h45m"]
        entries2 << [disk:"4GB"]
        entries2 << [maxForks:3]
        def config2 = new TaskConfig(entries2)
        def task2 = [name: 'Hello 2', processor: processor, config: config2, workDir: workDir2] as TaskRun
        def bean2 = new TaskBean(workDir: workDir2)
        def strategy2 = new WrFileCopyStrategy(bean2, null)

        def HTTPS_CONN = Mock(HttpsURLConnection)

        when:
        def jobs = client.add([
            ["exe one", task1, strategy1],
            ["exe two", task2, strategy2],
        ])

        then:
        1 * client.postJson('/jobs/', [
            [cmd:'exe one', cwd:'/my/work/dir', cwd_matters:true, rep_grp:'[nextflow] myprocess', req_grp:'[nextflow] myprocess', limit_grps:[], override:0, retries:0, rerun:true, cpus:1, memory:'', time:'', disk:0, mounts:[], on_exit:[[cleanup:true]]],
            [cmd:'exe two', cwd:'/tmp', cwd_matters:false, rep_grp:'[nextflow] myprocess', req_grp:'[nextflow] myprocess', limit_grps:["myprocess:3"], override:2, retries:0, rerun:true, cpus:2, memory:'500M', time:'2h', disk:4, mounts:[[Mount:'.mnt/bucket/work', Targets:[[Path:'bucket/work', Write:true]]]], on_exit:[['cleanup':true]]]
        ])
        1 * client.openConnection('https://localhost:12345/rest/v1/jobs/') >> HTTPS_CONN
        1 * client.authenticatedConnection('/jobs/')
        1 * HTTPS_CONN.setRequestProperty("Content-Type", "application/json") >> null
        1 * HTTPS_CONN.setRequestProperty("Authorization", "Bearer token") >> null
        1 * client.postConnection('/jobs/')
        1 * HTTPS_CONN.setRequestMethod('POST') >> null
        1 * HTTPS_CONN.setDoOutput(true) >> null
        1 * HTTPS_CONN.getOutputStream() >> new OutputStream() {
            @Override
            void write(int i) throws IOException { }
        }
        1 * HTTPS_CONN.getInputStream() >> { new ByteArrayInputStream('[{"k1":"v1"},{"k2":"v2"}]'.bytes) }
        1 * HTTPS_CONN.getResponseCode() >> 201
        jobs.size() == 2
    }

    def 'should get status' () {
        given:
        def client = Spy(createClient())
        def HTTPS_CONN = Mock(HttpsURLConnection)

        when:
        def jobs = client.status("a,b")

        then:
        1 * client.openConnection('https://localhost:12345/rest/v1/jobs/a,b') >> HTTPS_CONN
        1 * client.authenticatedConnection('/jobs/a,b')
        1 * HTTPS_CONN.setRequestProperty("Content-Type", "application/json") >> null
        1 * HTTPS_CONN.setRequestProperty("Authorization", "Bearer token") >> null
        1 * HTTPS_CONN.setRequestMethod('GET') >> null
        1 * HTTPS_CONN.getInputStream() >> { new ByteArrayInputStream('[{"k1":"v1"},{"k2":"v2"}]'.bytes) }
        1 * HTTPS_CONN.getResponseCode() >> 200
        jobs.size() == 2
    }

    def 'should cancel tasks in various states' () {
        given:
        def client = Spy(createClient())
        def HTTPS_CONN = Mock(HttpsURLConnection)

        when:
        client.cancel("a")

        then:
        1 * client.status("a") >> [[State:'buried']]
        1 * client.openConnection('https://localhost:12345/rest/v1/jobs/a?state=deletable') >> HTTPS_CONN
        1 * client.authenticatedConnection('/jobs/a?state=deletable')
        1 * HTTPS_CONN.setRequestProperty("Content-Type", "application/json") >> null
        1 * HTTPS_CONN.setRequestProperty("Authorization", "Bearer token") >> null
        1 * HTTPS_CONN.setRequestMethod('DELETE') >> null
        1 * HTTPS_CONN.getResponseCode() >> 200

        when:
        client.cancel("b")

        then:
        1 * client.status("b") >> [[State:'running']]
        1 * client.openConnection('https://localhost:12345/rest/v1/jobs/b?state=running') >> HTTPS_CONN
        1 * client.authenticatedConnection('/jobs/b?state=running')
        2 * HTTPS_CONN.setRequestProperty("Content-Type", "application/json") >> null
        2 * HTTPS_CONN.setRequestProperty("Authorization", "Bearer token") >> null
        2 * HTTPS_CONN.setRequestMethod('DELETE') >> null
        1 * HTTPS_CONN.getResponseCode() >> 202
        1 * client.status("b") >> [[State:'buried']]
        1 * client.openConnection('https://localhost:12345/rest/v1/jobs/b?state=deletable') >> HTTPS_CONN
        1 * client.authenticatedConnection('/jobs/b?state=deletable')
        1 * HTTPS_CONN.getResponseCode() >> 200


        when:
        client.cancel("c")

        then:
        1 * client.status("c") >> [[State:'lost']]
        1 * client.openConnection('https://localhost:12345/rest/v1/jobs/c?state=lost') >> HTTPS_CONN
        1 * client.authenticatedConnection('/jobs/c?state=lost')
        2 * HTTPS_CONN.setRequestProperty("Content-Type", "application/json") >> null
        2 * HTTPS_CONN.setRequestProperty("Authorization", "Bearer token") >> null
        2 * HTTPS_CONN.setRequestMethod('DELETE') >> null
        1 * HTTPS_CONN.getResponseCode() >> 202
        1 * client.status("c") >> [[State:'buried']]
        1 * client.openConnection('https://localhost:12345/rest/v1/jobs/c?state=deletable') >> HTTPS_CONN
        1 * client.authenticatedConnection('/jobs/c?state=deletable')
        1 * HTTPS_CONN.getResponseCode() >> 200
    }

}
