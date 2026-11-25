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

package nextflow.k8s.client

import java.nio.file.Files

import nextflow.exception.K8sOutOfCpuException
import nextflow.exception.K8sOutOfMemoryException

import javax.net.ssl.HttpsURLConnection

import nextflow.exception.NodeTerminationException
import nextflow.exception.ProcessFailedException
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sClientTest extends Specification {

    def 'should create a request' () {

        given:
        final TOKEN = '8d09d0ds'
        final client = Spy(K8sClient)

        def HTTPS_CONN = Mock(HttpsURLConnection)
        def HTTP_CONN = Mock(HttpURLConnection)

        when:
        client.config.server = 'host.com:443'
        client.config.token = TOKEN
        def resp = client.makeRequest('GET', '/foo/bar')
        then:
        1 * client.createConnection0("https://host.com:443/foo/bar") >> HTTPS_CONN
        1 * client.setupHttpsConn(HTTPS_CONN) >> { }
        1 * HTTPS_CONN.setRequestMethod('GET') >> { }
        1 * HTTPS_CONN.setRequestProperty("Authorization", "Bearer $TOKEN")
        1 * HTTPS_CONN.setRequestProperty("Content-Type", "application/json")
        1 * HTTPS_CONN.getResponseCode() >> 200
        1 * HTTPS_CONN.getInputStream() >> { new ByteArrayInputStream('{"field_x":"hello"}'.bytes) }
        resp instanceof K8sResponseApi
        resp.text == '{"field_x":"hello"}'

        when:
        client.config.server = 'http://my-server.com'
        client.config.token = TOKEN
        client.makeRequest('POST', '/foo/bar')
        then:
        1 * client.createConnection0("http://my-server.com/foo/bar") >> HTTP_CONN
        0 * client.setupHttpsConn(_) >> { }
        1 * HTTP_CONN.setRequestMethod('POST') >> { }
        1 * HTTP_CONN.getResponseCode() >> 401
        1 * HTTP_CONN.getErrorStream() >> { new ByteArrayInputStream('{"field_x":"oops.."}'.bytes) }
        def e = thrown(K8sResponseException)
        e.response.field_x == 'oops..'

    }

    def 'should make a get request' () {

        given:
        def client = Spy(K8sClient)
        when:
        client.get('/foo/bar')
        then:
        1 * client.makeRequest('GET', '/foo/bar') >> null
    }

    def 'should make a post request' () {

        given:
        def client = Spy(K8sClient)
        when:
        client.post('/foo/bar', '{ the: body }')
        then:
        1 * client.makeRequest('POST', '/foo/bar', '{ the: body }') >> null
    }

    def 'should make a delete request' () {

        given:
        def client = Spy(K8sClient)
        when:
        client.delete('/foo/bar', '{ the: body }')
        then:
        1 * client.makeRequest('DELETE', '/foo/bar', '{ the: body }') >> null
    }

    def 'should delete a pod' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"field":"OK"}'

        when:
        result = client.podDelete('foo')
        then:
        1 * client.delete("/api/v1/namespaces/default/pods/foo",null) >> RESP
        result.field == "OK"

        when:
        client.config.namespace = 'bar'
        result = client.podDelete('foo')
        then:
        1 * client.delete("/api/v1/namespaces/bar/pods/foo",null) >> RESP
        result.field == "OK"

    }

    def 'should list pods' () {
        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"hello"}'

        when:
        result = client.podList(true)
        then:
        1 * client.get('/api/v1/pods') >> RESP
        result.response == 'hello'

        when:
        result = client.podList()
        then:
        1 * client.get('/api/v1/namespaces/default/pods') >> RESP
        result.response == 'hello'

        when:
        client.config.namespace = 'foo'
        result = client.podList()
        then:
        1 * client.get('/api/v1/namespaces/foo/pods') >> RESP
        result.response == 'hello'
    }

    def 'should list secrets' () {
        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"hello"}'

        when:
        result = client.secretesList()
        then:
        1 * client.get("/api/v1/namespaces/default/secrets") >> RESP
        result.response == 'hello'

        when:
        client.config.namespace = 'pippo'
        result = client.podList()
        then:
        1 * client.get('/api/v1/namespaces/pippo/pods') >> RESP
        result.response == 'hello'

    }

    def 'should describe secret' () {
        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"hello"}'

        when:
        result = client.secretDescribe('my-secret')
        then:
        1 * client.get("/api/v1/namespaces/default/secrets/my-secret") >> RESP
        result.response == 'hello'

        when:
        client.config.namespace = 'paperino'
        result = client.secretDescribe('pluto')
        then:
        1 * client.get('/api/v1/namespaces/paperino/secrets/pluto') >> RESP
        result.response == 'hello'
    }

    def 'should log pod' () {

        given:
        InputStream result
        def client = Spy(K8sClient)
        def STREAM = Mock(InputStream)
        def RESP = Mock(K8sResponseApi)
        RESP.getStream() >> STREAM

        when:
        result = client.podLog('pod-123')
        then:
        1 * client.get('/api/v1/namespaces/default/pods/pod-123/log') >> RESP
        result == STREAM

        when:
        result = client.podLog('pod-123', follow: true)
        then:
        1 * client.get('/api/v1/namespaces/default/pods/pod-123/log?follow=true') >> RESP
        result == STREAM

        when:
        result = client.podLog('pod-123', follow: true, foo:1, bar: 'x')
        then:
        1 * client.get('/api/v1/namespaces/default/pods/pod-123/log?follow=true&foo=1&bar=x') >> RESP
        result == STREAM

    }

    def 'should create config' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"done"}'
        def CONFIG = [foo: 'hello', bar:'world']
        def JSON1 = '{"apiVersion":"v1","kind":"ConfigMap","metadata":{"name":"foo","namespace":"default"},"data":{"foo":"hello","bar":"world"}}'
        def JSON2 = '{"apiVersion":"v1","kind":"ConfigMap","metadata":{"name":"foo","namespace":"bar"},"data":{"foo":"hello","bar":"world"}}'

        when:
        result = client.configCreate('foo' , CONFIG)
        then:
        1 * client.post("/api/v1/namespaces/default/configmaps", JSON1) >> RESP
        result.response == 'done'

        when:
        client.config.namespace = 'bar'
        result = client.configCreate('foo' , CONFIG)
        then:
        1 * client.post("/api/v1/namespaces/bar/configmaps", JSON2) >> RESP
        result.response == 'done'

    }

    def 'should delete config' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"done"}'

        when:
        result = client.configDelete('foo')
        then:
        1 * client.delete("/api/v1/namespaces/default/configmaps/foo") >> RESP
        result.response == 'done'

        when:
        client.config.namespace = 'ns-1'
        result = client.configDelete('foo')
        then:
        1 * client.delete("/api/v1/namespaces/ns-1/configmaps/foo") >> RESP
        result.response == 'done'

    }

    def 'should delete all configs' () {

        given:
        K8sResponseJson result
        def client = Spy(K8sClient)
        def RESP = Mock(K8sResponseApi)
        RESP.getText() >> '{"response":"done"}'

        when:
        result = client.configDeleteAll()
        then:
        1 * client.delete("/api/v1/namespaces/default/configmaps") >> RESP
        result.response == 'done'

        when:
        client.config.namespace = 'ns-1'
        result = client.configDeleteAll()
        then:
        1 * client.delete("/api/v1/namespaces/ns-1/configmaps") >> RESP
        result.response == 'done'

    }

    def 'should get a pod state' () {

        given:
        def JSON = '''
          {
             "kind": "Pod",
             "apiVersion": "v1",
             "metadata": {
                 "name": "pod-xyz",
                 "namespace": "default",
                 "selfLink": "/api/v1/namespaces/default/pods/pod-xyz/status",
                 "uid": "33390e2c-f84b-11e7-a89d-025000000001",
                 "resourceVersion": "119932",
                 "creationTimestamp": "2018-01-13T10:19:12Z",
                 "labels": {
                     "app": "nextflow"
                 }
             },


             "status": {
                 "phase": "Succeeded",
                 "conditions": [
                     {
                         "type": "Initialized",
                         "status": "True",
                         "lastProbeTime": null,
                         "lastTransitionTime": "2018-01-13T10:19:12Z",
                         "reason": "PodCompleted"
                     },
                     {
                         "type": "Ready",
                         "status": "False",
                         "lastProbeTime": null,
                         "lastTransitionTime": "2018-01-13T10:19:37Z",
                         "reason": "PodCompleted"
                     },
                     {
                         "type": "PodScheduled",
                         "status": "True",
                         "lastProbeTime": null,
                         "lastTransitionTime": "2018-01-13T10:19:12Z"
                     }
                 ],
                 "hostIP": "192.168.65.3",
                 "podIP": "10.1.0.25",
                 "startTime": "2018-01-13T10:19:12Z",
                 "containerStatuses": [
                     {
                         "name": "pod-xyz",
                         "state": {
                             "terminated": {
                                 "exitCode": 0,
                                 "reason": "Completed",
                                 "startedAt": "2018-01-13T10:19:16Z",
                                 "finishedAt": "2018-01-13T10:19:36Z",
                                 "containerID": "docker://90d447d4a2518642b12c8979474aa2bbb5fe8e96ed9e5caf3c979bb3b751519a"
                             }
                         },
                         "lastState": {

                         },
                         "ready": false,
                         "restartCount": 0,
                         "image": "debian:latest",
                         "imageID": "docker-pullable://debian@sha256:0a5fcee6f52d5170f557ee2447d7a10a5bdcf715dd7f0250be0b678c556a501b",
                         "containerID": "docker://90d447d4a2518642b12c8979474aa2bbb5fe8e96ed9e5caf3c979bb3b751519a"
                     }
                 ],
                 "qosClass": "BestEffort"
             }
         }
        '''

        def client = Spy(K8sClient)
        final POD_NAME = 'pod-xyz'

        when:
        def result = client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)

        result == [terminated: [exitCode:0,
                                reason: 'Completed',
                                startedAt: "2018-01-13T10:19:16Z",
                                finishedAt: "2018-01-13T10:19:36Z",
                                containerID: "docker://90d447d4a2518642b12c8979474aa2bbb5fe8e96ed9e5caf3c979bb3b751519a"]]


    }

    def 'should get a pod its node' () {

        given:
        def JSON = '''
          {
             "kind": "Pod",
             "apiVersion": "v1",
             "metadata": {
                 "name": "pod-xyz",
                 "namespace": "default",
                 "selfLink": "/api/v1/namespaces/default/pods/pod-xyz/status",
                 "uid": "33390e2c-f84b-11e7-a89d-025000000001",
                 "resourceVersion": "119932",
                 "creationTimestamp": "2018-01-13T10:19:12Z",
                 "labels": {
                     "app": "nextflow"
                 }
             },
             "spec": {
                  "restartPolicy": "Never",
                  "terminationGracePeriodSeconds": 30,
                  "dnsPolicy": "ClusterFirst",
                  "serviceAccountName": "default",
                  "serviceAccount": "default",
                  "nodeName": "gke-bioinformatics-s-pipeline-pool-sm-bbac2e1c-k1tw",
                  "priority": 0,
                  "enableServiceLinks": true,
                  "preemptionPolicy": "PreemptLowerPriority"
              }
         }
        '''

        def client = Spy(K8sClient)
        final POD_NAME = 'pod-xyz'

        when:
        def result = client.getNodeOfPod(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)

        result == "gke-bioinformatics-s-pipeline-pool-sm-bbac2e1c-k1tw"

    }

    def 'should return undetermined status' () {
        given:
        def JSON = '''
             {
                 "kind": "Pod",
                 "apiVersion": "v1",
                 "metadata": {
                     "name": "nf-eb853c8010b8e173b23d8d15489d1a31",
                     "namespace": "default",
                     "selfLink": "/api/v1/namespaces/default/pods/nf-eb853c8010b8e173b23d8d15489d1a31/status",
                     "uid": "5ccdeb23-4f69-11e8-89b1-fa163e31bb09",
                     "resourceVersion": "2847182",
                     "creationTimestamp": "2018-05-04T07:04:18Z",
                     "labels": {
                         "app": "nextflow",
                         "processName": "markDuplicates",
                         "runName": "grave-jones",
                         "sessionId": "uuid-f51cd941-c21c-447b-86ca-eaebafa5ad9b",
                         "taskName": "markDuplicates_22028_2_118_1AlignedByCoord.out"
                     }
                 },

                 "status": {
                     "phase": "Pending",
                     "conditions": [
                         {
                             "type": "PodScheduled",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-04T07:04:18Z"
                         }
                     ],
                     "qosClass": "Guaranteed"
                 }
             }
        '''

        def client = Spy(K8sClient)
        final POD_NAME = 'nf-eb853c8010b8e173b23d8d15489d1a31'

        when:
        def result = client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)

        result == [:]
    }

    def 'should return undetermined status when status conditions are missing' () {
        given:
        def JSON = '''
            {
                "kind": "Pod",
                "apiVersion": "v1",
                "metadata": {
                    "name": "nf-89b34d7c2d1dc11daad72b1fcf7e0540",
                    "namespace": "default",
                    "selfLink": "/api/v1/namespaces/default/pods/nf-89b34d7c2d1dc11daad72b1fcf7e0540/status",
                    "uid": "02c3a34d-7720-11e9-8e1a-0a8b849038f8",
                    "resourceVersion": "21753960",
                    "creationTimestamp": "2019-05-15T14:44:58Z",
                    "labels": {
                        "app": "nextflow",
                        "processName": "split_input_file",
                        "runName": "hopeful_lalande",
                        "sessionId": "uuid-6b10c771-1a70-4af0-92be-7c0280e0e17f",
                        "taskName": "split_input_file_variants.gord_chr11_72000001-78000000",
                    }
                },

                "status": {
                    "phase": "Pending",
                    "qosClass": "Guaranteed"
                }
            }
        '''

        def client = Spy(K8sClient)
        final POD_NAME = 'nf-89b34d7c2d1dc11daad72b1fcf7e0540'

        when:
        def result = client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)

        result == [:]
    }

    def 'should return a process execution on pod not found' () {
        given:
        def JSON = '''
               {
                  "kind": "Status",
                  "apiVersion": "v1",
                  "metadata": {

                  },
                  "status": "Failure",
                  "message": "pods \\"nf-7cee928c1dd05b39cd50ab79a3e742f9\\" not found",
                  "reason": "NotFound",
                  "details": {
                      "name": "nf-7cee928c1dd05b39cd50ab79a3e742f9",
                      "kind": "pods"
                  },
                  "code": 404
              }
        '''

        def client = Spy(K8sClient)
        final POD_NAME = 'pod-xyz'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> { throw new K8sResponseException("Request GET /api/v1/namespaces/xyz/pods/nf-xyz/status returned an error code=404", new ByteArrayInputStream(JSON.bytes)) }

        and:
        thrown(NodeTerminationException)
    }

    def 'should fail to get pod state' () {

        given:
        def client = Spy(K8sClient)
        final POD_NAME = 'pod-xyz'
        final STATE = [foo:1, bar:2]

        when:
        def e
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson([:])
        e = thrown(K8sResponseException)
        e.message.startsWith('K8s undetermined status conditions for pod')

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson([status:[containerStatuses: []]])
        e = thrown(K8sResponseException)
        e.message.startsWith('K8s undetermined status conditions for pod')

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson([status:[containerStatuses: [ [name: 'foo'] ]]])
        e = thrown(K8sResponseException)
        e.message.startsWith('K8s invalid status for pod: pod-xyz (unexpected container name: foo)')

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson([status:[containerStatuses: [ [name: POD_NAME] ]]])
        e = thrown(K8sResponseException)
        e.message.startsWith('K8s invalid status for pod: pod-xyz (missing state object)')

        when:
        def result = client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson([status:[containerStatuses: [ [name: POD_NAME, state: STATE] ]]])
        result == STATE
    }

    def 'client should throw an exception when container status returns ErrImagePull' () {
        given:
        def JSON = '''
             {
                 "kind": "Pod",
                 "apiVersion": "v1",
                 "metadata": {
                     "name": "nf-c6e49ee9ebc79486f774a47924a743d7",
                     "namespace": "default",
                     "selfLink": "/api/v1/namespaces/default/pods/nf-c6e49ee9ebc79486f774a47924a743d7/status",
                     "uid": "18954b48-5811-11e8-8e71-025000000001",
                     "resourceVersion": "3615842",
                     "creationTimestamp": "2018-05-15T07:25:08Z",
                     "labels": {
                         "app": "nextflow",
                         "processName": "sayHello",
                         "runName": "lethal-jones",
                         "sessionId": "uuid-b6243a3d-5c07-44f4-97eb-9de499747800",
                         "taskName": "sayHello_2"
                     }
                 },
                 "spec": {

                 },
                 "status": {
                     "phase": "Pending",
                     "conditions": [
                         {
                             "type": "Initialized",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z"
                         },
                         {
                             "type": "Ready",
                             "status": "False",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z",
                             "reason": "ContainersNotReady",
                             "message": "containers with unready status: [nf-c6e49ee9ebc79486f774a47924a743d7]"
                         },
                         {
                             "type": "PodScheduled",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z"
                         }
                     ],
                     "hostIP": "192.168.65.3",
                     "podIP": "10.1.3.173",
                     "startTime": "2018-05-15T07:25:08Z",
                     "containerStatuses": [
                         {
                             "name": "nf-c6e49ee9ebc79486f774a47924a743d7",
                             "state": {
                                 "waiting": {
                                     "reason": "ErrImagePull",
                                     "message": "rpc error: code = Unknown desc = Error response from daemon: pull access denied for nextflow/foo, repository does not exist or may require 'docker login'"
                                 }
                             },
                             "lastState": {

                             },
                             "ready": false,
                             "restartCount": 0,
                             "image": "nextflow/foo",
                             "imageID": ""
                         }
                     ],
                     "qosClass": "BestEffort"
                 }
             }
'''

        def client = Spy(K8sClient)
        final POD_NAME = 'nf-c6e49ee9ebc79486f774a47924a743d7'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)
        def e = thrown(PodUnschedulableException)
        e.message == "K8s pod image cannot be pulled -- rpc error: code = Unknown desc = Error response from daemon: pull access denied for nextflow/foo, repository does not exist or may require 'docker login'"
    }

    def 'client should throw an exception when k8s is out of cpu' () {
        given:
        def JSON = '''
             {
                "kind": "Pod",
                "apiVersion": "v1",
                "metadata": {
                    "name": "nf-3b344812fe0aeb9554424bcf6caa7ffb",
                    "namespace": "default",
                    "uid": "0dd3c071-82a3-4b20-bdd4-0d34ff3d90bd",
                    "resourceVersion": "55320",
                    "creationTimestamp": "2022-09-23T13:43:48Z",
                    "labels": {
                        "app": "nextflow",
                        "processName": "combineFiles",
                        "runName": "insane-kare",
                        "sessionId": "uuid-85951b91-b5bc-4566-8938-fef4256c21c8",
                        "taskName": "combineFiles_1"
                    },
                },
                "spec": {
                },
                "status": {
                    "phase": "Failed",
                    "message": "Pod Node didn't have enough resource: cpu, requested: 4000, used: 2100, capacity: 6000",
                    "reason": "OutOfcpu",
                    "startTime": "2022-09-23T13:43:48Z"
                }
            }
        '''

        def client = Spy(K8sClient)
        final POD_NAME = 'nf-3b344812fe0aeb9554424bcf6caa7ffb'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)
        def e = thrown(K8sOutOfCpuException)
        e.message == "K8s pod 'nf-3b344812fe0aeb9554424bcf6caa7ffb' execution failed - reason: OutOfcpu - message: Pod Node didn't have enough resource: cpu, requested: 4000, used: 2100, capacity: 6000"
    }

    def 'client should throw an exception when k8s is out of memory' () {
        given:
        def JSON = '''
             {
                "kind": "Pod",
                "apiVersion": "v1",
                "metadata": {
                    "name": "nf-3b344812fe0aeb9554424bcf6caa7ffb",
                    "namespace": "default",
                    "uid": "0dd3c071-82a3-4b20-bdd4-0d34ff3d90bd",
                    "resourceVersion": "55320",
                    "creationTimestamp": "2022-09-23T13:43:48Z",
                    "labels": {
                        "app": "nextflow",
                        "processName": "combineFiles",
                        "runName": "insane-kare",
                        "sessionId": "uuid-85951b91-b5bc-4566-8938-fef4256c21c8",
                        "taskName": "combineFiles_2"
                    },
                },
                "spec": {
                },
                "status": {
                    "phase": "Failed",
                    "message": "Pod Node didn't have enough resource: memory, requested: 16106127360, used: 16158556160, capacity: 16778358784",
                    "reason": "OutOfmemory",
                    "startTime": "2022-09-23T13:42:59Z"
                }
            }
        '''

        def client = Spy(K8sClient)
        final POD_NAME = 'nf-3b344812fe0aeb9554424bcf6caa7ffb'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)
        def e = thrown(K8sOutOfMemoryException)
        e.message == "K8s pod 'nf-3b344812fe0aeb9554424bcf6caa7ffb' execution failed - reason: OutOfmemory - message: Pod Node didn't have enough resource: memory, requested: 16106127360, used: 16158556160, capacity: 16778358784"
    }

    def 'client should throw an exception when container status returns ImagePullBackOff' () {
        def JSON = '''
            {
                 "kind": "Pod",
                 "apiVersion": "v1",
                 "metadata": {
                     "name": "nf-c6e49ee9ebc79486f774a47924a743d7",
                     "namespace": "default",
                     "selfLink": "/api/v1/namespaces/default/pods/nf-c6e49ee9ebc79486f774a47924a743d7/status",
                     "uid": "18954b48-5811-11e8-8e71-025000000001",
                     "resourceVersion": "3615866",
                     "creationTimestamp": "2018-05-15T07:25:08Z",
                     "labels": {
                         "app": "nextflow",
                         "processName": "sayHello",
                         "runName": "lethal-jones",
                         "sessionId": "uuid-b6243a3d-5c07-44f4-97eb-9de499747800",
                         "taskName": "sayHello_2"
                     }
                 },
                 "spec": { },
                 "status": {
                     "phase": "Pending",
                     "conditions": [
                         {
                             "type": "Initialized",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z"
                         },
                         {
                             "type": "Ready",
                             "status": "False",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z",
                             "reason": "ContainersNotReady",
                             "message": "containers with unready status: [nf-c6e49ee9ebc79486f774a47924a743d7]"
                         },
                         {
                             "type": "PodScheduled",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-05-15T07:25:08Z"
                         }
                     ],
                     "hostIP": "192.168.65.3",
                     "podIP": "10.1.3.173",
                     "startTime": "2018-05-15T07:25:08Z",
                     "containerStatuses": [
                         {
                             "name": "nf-c6e49ee9ebc79486f774a47924a743d7",
                             "state": {
                                 "waiting": {
                                     "reason": "ImagePullBackOff",
                                     "message": "Back-off pulling image \\"nextflow/foo\\""
                                 }
                             },
                             "lastState": {

                             },
                             "ready": false,
                             "restartCount": 0,
                             "image": "nextflow/foo",
                             "imageID": ""
                         }
                     ],
                     "qosClass": "BestEffort"
                 }
             }
'''
        def client = Spy(K8sClient)
        final POD_NAME = 'nf-c6e49ee9ebc79486f774a47924a743d7'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)
        def e = thrown(PodUnschedulableException)
        e.message == "K8s pod image cannot be pulled -- Back-off pulling image \"nextflow/foo\""

    }

    def 'client should throw process exception on failed state' () {
        def JSON = '''
             {
              "kind": "Pod",
              "apiVersion": "v1",
              "metadata": {
                  "name": "nf-f34e124e1471736f27c8ef1aa52d02be",
                  "namespace": "tower-nf",
                  "uid": "8ddb536f-da4c-4a99-a5c5-3b5b53df6411",
                  "resourceVersion": "31269576",
                  "creationTimestamp": "2021-08-09T17:58:42Z",
                  "labels": {
                      "app": "nextflow",
                      "processName": "bulk_rnaseq_trim_galore",
                      "runName": "insane_agnesi",
                      "sessionId": "uuid-35a8fe1a-2622-4582-a1b2-1cc21bd85b57",
                      "taskName": "bulk_rnaseq_trim_galore_d0_tdTOM-0701"
                  }
              },
              "spec": {
                  "restartPolicy": "Never",
                  "terminationGracePeriodSeconds": 30,
                  "dnsPolicy": "ClusterFirst",
                  "serviceAccountName": "default",
                  "serviceAccount": "default",
                  "nodeName": "gke-bioinformatics-s-pipeline-pool-sm-bbac2e1c-k1tw",
                  "priority": 0,
                  "enableServiceLinks": true,
                  "preemptionPolicy": "PreemptLowerPriority"
              },
              "status": {
                  "phase": "Failed",
                  "message": "Node is shutting, evicting pods",
                  "reason": "Shutdown",
                  "startTime": "2021-08-09T17:59:31Z"
              }
            }
'''
        def client = Spy(K8sClient)
        final POD_NAME = 'nf-xyz'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)
        and:
        def e = thrown(NodeTerminationException)
        e.message == "K8s pod 'nf-xyz' execution failed - reason: Shutdown - message: Node is shutting, evicting pods"

    }

    def 'client should fail when config fail' () {
        given:
        def JSON = '''
             {
                 "kind": "Pod",
                 "apiVersion": "v1",
                 "metadata": {
                     "name": "angry-blackwell",
                     "namespace": "default",
                     "selfLink": "/api/v1/namespaces/default/pods/angry-blackwell/status",
                     "uid": "83a35c1e-73b6-11e8-8259-025000000001",
                     "resourceVersion": "465382",
                     "creationTimestamp": "2018-06-19T11:47:16Z",
                     "labels": {
                         "app": "nextflow",
                         "runName": "angry-blackwell"
                     }
                 },
                 "spec": {

                 },
                 "status": {
                     "phase": "Pending",
                     "conditions": [
                         {
                             "type": "Initialized",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-06-19T11:47:16Z"
                         },
                         {
                             "type": "Ready",
                             "status": "False",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-06-19T11:47:16Z",
                             "reason": "ContainersNotReady",
                             "message": "containers with unready status: [angry-blackwell]"
                         },
                         {
                             "type": "PodScheduled",
                             "status": "True",
                             "lastProbeTime": null,
                             "lastTransitionTime": "2018-06-19T11:47:16Z"
                         }
                     ],
                     "hostIP": "192.168.65.3",
                     "podIP": "10.1.4.20",
                     "startTime": "2018-06-19T11:47:16Z",
                     "containerStatuses": [
                         {
                             "name": "angry-blackwell",
                             "state": {
                                 "waiting": {
                                     "reason": "CreateContainerConfigError",
                                     "message": "secrets \\"my-env\\" not found"
                                 }
                             },
                             "lastState": {

                             },
                             "ready": false,
                             "restartCount": 0,
                             "image": "nextflow/nextflow:0.30.2",
                             "imageID": ""
                         }
                     ],
                     "qosClass": "BestEffort"
                 }
             }
            '''

        def client = Spy(K8sClient)
        final POD_NAME = 'angry-blackwell'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)
        def e = thrown(PodUnschedulableException)
        e.message == 'K8s pod configuration failed -- secrets "my-env" not found'

    }

    def 'client should throw process exception when initContainer fails' () {
        // Addresses bug 2428, which is when an initContainer exists and fails
        // leaving the containerStatuses in "PodInitializing", but the pod
        // phase is "Failed".
        def JSON = '''
            {
                "apiVersion": "v1",
                "kind": "Pod",
                "metadata": {
                    "creationTimestamp": "2021-11-03T21:37:45Z",
                    "labels": {
                        "app": "nextflow",
                        "magic-version": "v7.5.3",
                        "processName": "some_special_magic",
                        "runName": "lonely_tuckerman",
                        "sessionId": "uuid-41ef6f42-bab5-438f-83e4-85e91a6386ea",
                        "taskName": "some_special_magic_1"
                    },
                    "name": "nf-45d62fd8462ba390bacfa20ad9065bfe",
                    "namespace": "default",
                    "resourceVersion": "90692",
                    "uid": "280bf29c-4c14-4451-a300-080760ff7a0e"
                },
                "status": {
                    "containerStatuses": [
                        {
                            "image": "docker/asdf/magictest:7.5.3",
                            "imageID": "",
                            "lastState": {},
                            "name": "nf-45d62fd8462ba390bacfa20ad9065bfe",
                            "ready": false,
                            "restartCount": 0,
                            "started": false,
                            "state": {
                                "waiting": {
                                    "reason": "PodInitializing"
                                }
                            }
                        }
                    ],
                    "hostIP": "10.12.205.49",
                    "initContainerStatuses": [
                        {
                            "containerID": "containerd://8ac143e547ef8c417f3655da31a81018e9d14db08c218e56ee3d181c70e0b755",
                            "image": "docker/asdf/magic-special-handling:0.4.0",
                            "imageID": "sha256:587f96e98fd0a2b17772e9a505522e4c0e1cf5194a4f08924c339848df42100c",
                            "lastState": {},
                            "name": "magic-init",
                            "ready": false,
                            "restartCount": 0,
                            "state": {
                                "terminated": {
                                    "containerID": "containerd://8ac143e547ef8c417f3655da31a81018e9d14db08c218e56ee3d181c70e0b755",
                                    "exitCode": 1,
                                    "finishedAt": "2021-11-03T21:37:51Z",
                                    "reason": "Error",
                                    "startedAt": "2021-11-03T21:37:46Z"
                                }
                            }
                        }
                    ],
                    "phase": "Failed",
                    "podIP": "10.42.0.67",
                    "podIPs": [
                        {
                            "ip": "10.42.0.67"
                        }
                    ],
                    "qosClass": "Burstable",
                    "startTime": "2021-11-03T21:37:45Z"
                }
            }
'''
        def client = Spy(K8sClient)
        final POD_NAME = 'nf-45d62fd8462ba390bacfa20ad9065bfe'

        when:
        client.podState(POD_NAME)
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(JSON)
        and:
        def e = thrown(PodUnschedulableException)
        e.message == "K8s pod in Failed state"
    }

    def 'should fallback to job status when pod is gone and not return hardcoded exit code' () {
        given:
        def JOB_STATUS_JSON = '''
        {
            "apiVersion": "batch/v1",
            "kind": "Job",
            "metadata": {
                "name": "test-job"
            },
            "status": {
                "succeeded": 1,
                "startTime": "2025-01-15T10:00:00Z",
                "completionTime": "2025-01-15T10:05:00Z",
                "conditions": [
                    {
                        "type": "Complete",
                        "status": "True",
                        "lastProbeTime": "2025-01-15T10:05:00Z",
                        "lastTransitionTime": "2025-01-15T10:05:00Z"
                    }
                ]
            }
        }
        '''
        def client = Spy(K8sClient)
        final JOB_NAME = 'test-job'

        when:
        def result = client.jobStateFallback0(JOB_NAME)

        then:
        1 * client.jobStatus(JOB_NAME) >> new K8sResponseJson(JOB_STATUS_JSON)

        and:
        result.terminated != null
        result.terminated.reason == 'Completed'
        result.terminated.startedAt == '2025-01-15T10:00:00Z'
        result.terminated.finishedAt == '2025-01-15T10:05:00Z'
        // The key assertion: exitCode should not be present (null) so fallback to .exitcode file works
        result.terminated.exitCode == null
        result.terminated.exitcode == null
    }

    /*
     * pod states captured on a real cluster (Kubernetes v1.37) while disrupting a pod with
     * the scheduler preemption, the eviction API, a NoExecute taint, the kubelet eviction, a
     * graceful node shutdown, a node lost without notice and a node deleted from the cluster.
     * Only the attributes read by the client are kept, and the names are replaced.
     * See the pull request discussion for the script used to reproduce them on a kind cluster
     */
    static final Map<String,String> POD_STATES = [
        'running': '''{"metadata":{"name":"nf-4f5d3c2b1a"},"status":{"phase":"Running","conditions":[{"type":"Ready","status":"True"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"running":{"startedAt":"2026-11-11T15:48:57Z"}}}]}}''',
        'succeeded': '''{"metadata":{"name":"nf-4f5d3c2b1a"},"status":{"phase":"Succeeded","conditions":[{"type":"Ready","status":"False","reason":"PodCompleted"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"terminated":{"exitCode":0,"finishedAt":"2026-11-11T15:44:55Z","reason":"Completed","startedAt":"2026-11-11T15:44:55Z"}}}]}}''',
        'exit-error': '''{"metadata":{"name":"nf-4f5d3c2b1a"},"status":{"phase":"Failed","conditions":[{"type":"Ready","status":"False","reason":"PodFailed"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"terminated":{"exitCode":3,"finishedAt":"2026-11-11T15:45:33Z","reason":"Error","startedAt":"2026-11-11T15:45:33Z"}}}]}}''',
        'oom': '''{"metadata":{"name":"nf-4f5d3c2b1a"},"status":{"phase":"Failed","conditions":[{"type":"Ready","status":"False","reason":"PodFailed"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"terminated":{"exitCode":137,"finishedAt":"2026-11-11T15:46:09Z","reason":"OOMKilled","startedAt":"2026-11-11T15:46:09Z"}}}]}}''',
        'kubelet-eviction': '''{"metadata":{"name":"nf-4f5d3c2b1a"},"status":{"phase":"Failed","reason":"Evicted","message":"Usage of EmptyDir volume \\"scratch\\" exceeds the limit \\"100Mi\\". ","conditions":[{"type":"Ready","status":"False","reason":"PodFailed"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"terminated":{"exitCode":137,"finishedAt":"2026-11-11T15:47:04Z","reason":"Error","startedAt":"2026-11-11T15:46:47Z"}}}]}}''',
        'api-eviction-terminating': '''{"metadata":{"name":"nf-4f5d3c2b1a","deletionTimestamp":"2026-11-11T15:48:39Z"},"status":{"phase":"Running","conditions":[{"type":"DisruptionTarget","status":"True","reason":"EvictionByEvictionAPI","message":"Eviction API: evicting"},{"type":"Ready","status":"True"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"running":{"startedAt":"2026-11-11T15:48:25Z"}}}]}}''',
        'preemption-terminating': '''{"metadata":{"name":"nf-4f5d3c2b1a","deletionTimestamp":"2026-11-11T15:49:11Z"},"status":{"phase":"Running","conditions":[{"type":"DisruptionTarget","status":"True","reason":"PreemptionByScheduler","message":"default-scheduler: preempting to accommodate a higher priority pod"},{"type":"Ready","status":"True"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"running":{"startedAt":"2026-11-11T15:48:57Z"}}}]}}''',
        'preemption-terminated': '''{"metadata":{"name":"nf-4f5d3c2b1a","deletionTimestamp":"2026-11-11T15:49:26Z"},"status":{"phase":"Failed","conditions":[{"type":"DisruptionTarget","status":"True","reason":"PreemptionByScheduler","message":"default-scheduler: preempting to accommodate a higher priority pod"},{"type":"Ready","status":"False","reason":"PodFailed"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"terminated":{"exitCode":137,"finishedAt":"2026-11-11T15:49:26Z","reason":"Error","startedAt":"2026-11-11T15:49:13Z"}}}]}}''',
        'taint-terminating': '''{"metadata":{"name":"nf-4f5d3c2b1a","deletionTimestamp":"2026-11-11T15:51:59Z"},"status":{"phase":"Running","conditions":[{"type":"DisruptionTarget","status":"True","reason":"DeletionByTaintManager","message":"Taint manager: deleting due to NoExecute taint"},{"type":"Ready","status":"True"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"running":{"startedAt":"2026-11-11T15:51:45Z"}}}]}}''',
        'taint-terminated': '''{"metadata":{"name":"nf-4f5d3c2b1a","deletionTimestamp":"2026-11-11T15:53:19Z"},"status":{"phase":"Failed","conditions":[{"type":"DisruptionTarget","status":"True","reason":"DeletionByTaintManager","message":"Taint manager: deleting due to NoExecute taint"},{"type":"Ready","status":"False","reason":"PodFailed"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"terminated":{"exitCode":137,"finishedAt":"2026-11-11T15:53:19Z","reason":"Error","startedAt":"2026-11-11T15:53:05Z"}}}]}}''',
        'node-shutdown': '''{"metadata":{"name":"nf-4f5d3c2b1a"},"status":{"phase":"Failed","reason":"Terminated","message":"Pod was terminated in response to imminent node shutdown.","conditions":[{"type":"DisruptionTarget","status":"True","reason":"TerminationByKubelet","message":"Pod was terminated in response to imminent node shutdown."},{"type":"Ready","status":"False","reason":"PodFailed"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"terminated":{"exitCode":137,"finishedAt":"2026-11-24T15:52:03Z","reason":"Error","startedAt":"2026-11-24T15:51:50Z"}}}]}}''',
        'node-lost-terminating': '''{"metadata":{"name":"nf-4f5d3c2b1a","deletionTimestamp":"2026-11-24T15:37:46Z"},"status":{"phase":"Running","conditions":[{"type":"Ready","status":"False"},{"type":"DisruptionTarget","status":"True","reason":"DeletionByTaintManager","message":"Taint manager: deleting due to NoExecute taint"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"running":{"startedAt":"2026-11-24T15:31:48Z"}}}]}}''',
        'node-delete-terminated': '''{"metadata":{"name":"nf-4f5d3c2b1a","deletionTimestamp":"2026-11-24T15:47:16Z"},"status":{"phase":"Failed","conditions":[{"type":"DisruptionTarget","status":"True","reason":"DeletionByPodGC","message":"PodGC: node no longer exists"},{"type":"Ready","status":"False","reason":"PodFailed"}],"containerStatuses":[{"name":"nf-4f5d3c2b1a","state":{"running":{"startedAt":"2026-11-24T15:46:43Z"}}}]}}''',
    ]

    def 'should classify the pod states captured on a real cluster' () {
        given:
        def client = Spy(K8sClient)
        final POD_NAME = 'nf-4f5d3c2b1a'

        when:
        def state = null
        def error = null
        try { state = client.podState(POD_NAME) }
        catch( Throwable t ) { error = t }
        then:
        1 * client.podStatus(POD_NAME) >> new K8sResponseJson(POD_STATES[CAPTURE])
        and:
        error?.getClass() == ERROR
        error?.message == MESSAGE
        error || state.keySet().first() == STATE

        where:
        CAPTURE                     | ERROR                     | STATE         | MESSAGE
        'running'                   | null                      | 'running'     | null
        'succeeded'                 | null                      | 'terminated'  | null
        // task failures: the container state is returned, the handler reports its reason (e.g. OOMKilled)
        'exit-error'                | null                      | 'terminated'  | null
        'oom'                       | null                      | 'terminated'  | null
        // infrastructure failures
        'kubelet-eviction'          | NodeTerminationException  | null          | "K8s pod 'nf-4f5d3c2b1a' execution failed - reason: Evicted - message: Usage of EmptyDir volume \"scratch\" exceeds the limit \"100Mi\". "
        'api-eviction-terminating'  | NodeTerminationException  | null          | "K8s pod 'nf-4f5d3c2b1a' execution failed - disruption: EvictionByEvictionAPI (Eviction API: evicting)"
        'preemption-terminating'    | NodeTerminationException  | null          | "K8s pod 'nf-4f5d3c2b1a' execution failed - disruption: PreemptionByScheduler (default-scheduler: preempting to accommodate a higher priority pod)"
        'preemption-terminated'     | NodeTerminationException  | null          | "K8s pod 'nf-4f5d3c2b1a' execution failed - disruption: PreemptionByScheduler (default-scheduler: preempting to accommodate a higher priority pod)"
        'taint-terminating'         | NodeTerminationException  | null          | "K8s pod 'nf-4f5d3c2b1a' execution failed - disruption: DeletionByTaintManager (Taint manager: deleting due to NoExecute taint)"
        'taint-terminated'          | NodeTerminationException  | null          | "K8s pod 'nf-4f5d3c2b1a' execution failed - disruption: DeletionByTaintManager (Taint manager: deleting due to NoExecute taint)"
        'node-shutdown'             | NodeTerminationException  | null          | "K8s pod 'nf-4f5d3c2b1a' execution failed - reason: Terminated - message: Pod was terminated in response to imminent node shutdown. - disruption: TerminationByKubelet (Pod was terminated in response to imminent node shutdown.)"
        'node-lost-terminating'     | NodeTerminationException  | null          | "K8s pod 'nf-4f5d3c2b1a' execution failed - disruption: DeletionByTaintManager (Taint manager: deleting due to NoExecute taint)"
        'node-delete-terminated'    | NodeTerminationException  | null          | "K8s pod 'nf-4f5d3c2b1a' execution failed - disruption: DeletionByPodGC (PodGC: node no longer exists)"
    }

    def 'should check infrastructure failure' () {
        given:
        def client = Spy(K8sClient)
        def status = [phase: PHASE, reason: REASON]
        if( DISRUPTION )
            status.conditions = [[type: 'DisruptionTarget', status: 'True', reason: DISRUPTION]]
        if( STATE != null )
            status.containerStatuses = [[name: 'nf-xyz', state: STATE]]
        def pod = [metadata: DELETING ? [deletionTimestamp: '2026-11-11T15:49:11Z'] : [:], status: status]

        when:
        def error = null
        try { client.checkInfraFailure('nf-xyz', pod) }
        catch( Throwable t ) { error = t }
        then:
        error?.getClass() == EXPECTED

        where:
        PHASE       | REASON        | DISRUPTION                | DELETING  | STATE                         | EXPECTED
        // being deleted because of a disruption, the container is still running during the grace period
        'Running'   | null          | 'PreemptionByScheduler'   | true      | [running: [:]]                | NodeTerminationException
        'Running'   | null          | 'EvictionByEvictionAPI'   | true      | [running: [:]]                | NodeTerminationException
        'Pending'   | null          | 'DeletionByTaintManager'  | true      | [waiting: [:]]                | NodeTerminationException
        // the condition is removed later when the pod is not deleted, therefore it's not enough
        'Running'   | null          | 'PreemptionByScheduler'   | false     | [running: [:]]                | null
        // deleted without a disruption (e.g. by nextflow itself)
        'Running'   | null          | null                      | true      | [running: [:]]                | null
        // terminated by a node disruption
        'Failed'    | null          | 'PreemptionByScheduler'   | true      | [terminated: [exitCode: 137]] | NodeTerminationException
        'Failed'    | null          | 'TerminationByKubelet'    | false     | [terminated: [exitCode: 143]] | NodeTerminationException
        // completed successfully before being disrupted
        'Succeeded' | null          | 'PreemptionByScheduler'   | true      | [terminated: [exitCode: 0]]   | null
        // pod failures reported by the pod reason
        'Failed'    | 'Evicted'     | null                      | false     | [terminated: [exitCode: 137]] | NodeTerminationException
        'Failed'    | 'Evicted'     | null                      | false     | null                          | NodeTerminationException
        'Failed'    | 'Preempting'  | null                      | false     | [terminated: [exitCode: 137]] | NodeTerminationException
        'Failed'    | 'Shutdown'    | null                      | false     | [terminated: [exitCode: 143]] | NodeTerminationException
        'Failed'    | 'Shutdown'    | null                      | false     | null                          | NodeTerminationException
        'Failed'    | 'OutOfcpu'    | null                      | false     | null                          | K8sOutOfCpuException
        'Failed'    | 'OutOfmemory' | null                      | false     | null                          | K8sOutOfMemoryException
        // task failures
        'Failed'    | null          | null                      | false     | [terminated: [exitCode: 1]]   | null
        'Failed'    | 'DeadlineExceeded' | null                 | false     | null                          | null
    }

    def 'should report node termination when the job pod is disrupted' () {
        given:
        def JOB_STATUS_JSON = """
        {
            "apiVersion": "batch/v1",
            "kind": "Job",
            "metadata": { "name": "test-job" },
            "status": $JOB_STATUS
        }
        """
        def client = Spy(K8sClient)
        final JOB_NAME = 'test-job'
        final POD_NAME = 'test-job-abc12'
        final ERROR = new NodeTerminationException("K8s pod '$POD_NAME' execution failed - disruption: PreemptionByScheduler")

        when:
        client.jobState(JOB_NAME)
        then:
        1 * client.findPodNameForJob(JOB_NAME) >> POD_NAME
        1 * client.podState(POD_NAME) >> { throw ERROR }
        1 * client.jobStatus(JOB_NAME) >> new K8sResponseJson(JOB_STATUS_JSON)
        and:
        def e = thrown(NodeTerminationException)
        e.is(ERROR)

        where:
        JOB_STATUS << [
            // the pod is terminating and the job did not fail yet
            '{"active": 1, "startTime": "2026-11-11T15:49:12Z"}',
            // the job failed because of the pod disruption
            '{"failed": 1, "startTime": "2026-11-11T15:49:12Z", "conditions": [{"type": "FailureTarget", "status": "True", "reason": "BackoffLimitExceeded", "message": "Job has reached the specified backoff limit"}]}',
        ]
    }

    def 'should report the job failure when the pod is gone (listed: #LISTED)' () {
        given:
        def JOB_STATUS_JSON = '''
        {
            "apiVersion": "batch/v1",
            "kind": "Job",
            "metadata": { "name": "test-job" },
            "status": {
                "failed": 1,
                "startTime": "2026-11-11T15:49:12Z",
                "conditions": [ {"type": "Failed", "status": "True", "reason": "BackoffLimitExceeded", "message": "Job has reached the specified backoff limit"} ]
            }
        }
        '''
        def client = Spy(K8sClient)
        final JOB_NAME = 'test-job'
        final POD_NAME = 'test-job-abc12'
        final MISSING = '{"kind":"Status","details":{"kind":"pods","name":"test-job-abc12"},"code":404}'

        when:
        client.jobState(JOB_NAME)
        then:
        1 * client.findPodNameForJob(JOB_NAME) >> (LISTED ? POD_NAME : null)
        (LISTED ? 1 : 0) * client.podStatus(POD_NAME) >> {
            throw new K8sResponseException('Pod not found', new ByteArrayInputStream(MISSING.bytes))
        }
        1 * client.jobStatus(JOB_NAME) >> new K8sResponseJson(JOB_STATUS_JSON)
        and:
        def e = thrown(ProcessFailedException)
        e.message == "K8s Job $JOB_NAME execution failed: Job has reached the specified backoff limit"

        where:
        LISTED << [false, true]
    }

    def 'should re-read token from disk and retry on 401 when tokenPath is set' () {

        given:
        def folder = Files.createTempDirectory('test')
        def tokenFile = folder.resolve('token')
        // file already holds the rotated token (kubelet has written it)
        tokenFile.text = 'fresh-token'

        final client = Spy(K8sClient)
        client.config.server = 'host.com:443'
        client.config.token = 'stale-token'
        client.config.tokenPath = tokenFile
        // shorten retry delay so the test is fast
        client.config.retryConfig.delay = nextflow.util.Duration.of('1ms')

        def CONN_401 = Mock(HttpsURLConnection)
        def CONN_200 = Mock(HttpsURLConnection)

        when:
        def resp = client.makeRequest('GET', '/api/v1/pods/foo/status')

        then:
        // first attempt: stale token, 401
        2 * client.createConnection0("https://host.com:443/api/v1/pods/foo/status") >>> [CONN_401, CONN_200]
        2 * client.setupHttpsConn(_) >> { }
        1 * CONN_401.setRequestProperty("Authorization", "Bearer stale-token")
        1 * CONN_401.setRequestMethod('GET')
        1 * CONN_401.setRequestProperty("Content-Type", "application/json")
        1 * CONN_401.getResponseCode() >> 401
        1 * CONN_401.getErrorStream() >> { new ByteArrayInputStream('{"kind":"Status","status":"Failure","message":"Unauthorized","code":401}'.bytes) }

        and:
        // second attempt: fresh token from disk, 200
        1 * CONN_200.setRequestProperty("Authorization", "Bearer fresh-token")
        1 * CONN_200.setRequestMethod('GET')
        1 * CONN_200.setRequestProperty("Content-Type", "application/json")
        1 * CONN_200.getResponseCode() >> 200
        1 * CONN_200.getInputStream() >> { new ByteArrayInputStream('{"ok":true}'.bytes) }

        and:
        resp.text == '{"ok":true}'
        client.config.token == 'fresh-token'

        cleanup:
        folder?.deleteDir()
    }

    def 'should not retry 401 when tokenPath is not set' () {

        given:
        final client = Spy(K8sClient)
        client.config.server = 'host.com:443'
        client.config.token = 'inline-token'
        // tokenPath is null — token can't be re-read from disk

        def CONN = Mock(HttpsURLConnection)

        when:
        client.makeRequest('GET', '/api/v1/pods/foo/status')

        then:
        // exactly one attempt — no retry because there's no file to re-read
        1 * client.createConnection0(_) >> CONN
        1 * client.setupHttpsConn(_) >> { }
        1 * CONN.setRequestMethod('GET')
        1 * CONN.setRequestProperty("Authorization", "Bearer inline-token")
        1 * CONN.setRequestProperty("Content-Type", "application/json")
        1 * CONN.getResponseCode() >> 401
        1 * CONN.getErrorStream() >> { new ByteArrayInputStream('{"kind":"Status","status":"Failure","message":"Unauthorized","code":401}'.bytes) }

        and:
        thrown(K8sResponseException)
    }
}
