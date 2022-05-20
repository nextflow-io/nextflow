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

package nextflow.k8s.model

import nextflow.executor.res.AcceleratorResource
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Lukas Hejtmanek <xhejtman@gmail.com>
 */
class JobSpecBuilderTest extends Specification {

    def setup() {
        PodSpecBuilder.VOLUMES.set(0)
    }


    def 'should create job spec' () {

        when:
        def specp = new PodSpecBuilder()
                .withContainerName('foo')
                .withImageName('busybox')
                .withCommand(['echo', 'hello'])
                .build()
        def specj = new JobSpecBuilder()
                    .withJobName('foo')
                    .withPodSpec(specp)
                    .build()

        then:
        specj ==  [ apiVersion: 'batch/v1',
                   kind: 'Job',
                   metadata: [name:'foo', namespace:'default'],
                   spec: [
                       backoffLimit: 0,                      
                       template: [
                           apiVersion: 'v1',
                           kind: 'Pod',
                           metadata: [:],
                           spec: [
                               restartPolicy:'Never',
                               containers:[
                                       [name:'foo',
                                        image:'busybox',
                                        command:['echo', 'hello'],
                                       ]
                               ]
                           ]
                       ]
                   ]
                 ]
    }

    def 'should create job spec with backofflimit' () {

        when:
        def specp = new PodSpecBuilder()
                .withContainerName('foo')
                .withImageName('busybox')
                .withCommand(['echo', 'hello'])
                .build()
        def specj = new JobSpecBuilder()
                    .withJobName('foo')
                    .withPodSpec(specp)
                    .withJobBackoffLimit(6)
                    .build()

        then:
        specj ==  [ apiVersion: 'batch/v1',
                   kind: 'Job',
                   metadata: [name:'foo', namespace:'default'],
                   spec: [
                       backoffLimit: 6,                      
                       template: [
                           apiVersion: 'v1',
                           kind: 'Pod',
                           metadata: [:],
                           spec: [
                               restartPolicy:'Never',
                               containers:[
                                       [name:'foo',
                                        image:'busybox',
                                        command:['echo', 'hello'],
                                       ]
                               ]
                           ]
                       ]
                   ]
                 ]
    }
}
