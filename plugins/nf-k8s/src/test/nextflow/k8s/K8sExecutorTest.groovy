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

package nextflow.k8s

import nextflow.k8s.client.K8sClient
import spock.lang.Specification

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sExecutorTest extends Specification {

    def 'should return the same k8s client instance on repeated calls' () {
        given:
        def CONFIG = new K8sConfig(
            client: [server: 'http://k8s-server'],
            namespace: 'test-ns',
            serviceAccount: 'test-sa'
        )
        and:
        def executor = Spy(K8sExecutor)
        executor.getK8sConfig() >> CONFIG
        executor.@client = new K8sClient(CONFIG.getClient())

        when:
        def client1 = executor.getClient()
        def client2 = executor.getClient()

        then:
        client1 instanceof K8sClient
        client2.is(client1)
    }

}
