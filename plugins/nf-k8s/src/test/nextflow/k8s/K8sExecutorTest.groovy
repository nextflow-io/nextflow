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

import java.util.concurrent.TimeUnit

import com.google.common.cache.CacheBuilder
import nextflow.k8s.client.ClientConfig
import nextflow.k8s.client.K8sClient
import spock.lang.Specification

/**
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class K8sExecutorTest extends Specification {

    def 'should cache k8s client and refresh after expiration' () {
        given:
        def CONFIG = new K8sConfig(
            client: [server: 'http://k8s-server'],
            namespace: 'test-ns',
            serviceAccount: 'test-sa',
            clientRefreshInterval: '100ms'
        )
        and:
        def executor = Spy(K8sExecutor)
        executor.getK8sConfig() >> CONFIG
        // use a short-lived cache for the test
        executor.@clientCache = CacheBuilder.newBuilder()
            .expireAfterWrite(100, TimeUnit.MILLISECONDS)
            .build()

        when: 'first call to getClient'
        def client1 = executor.getClient()
        then: 'a new K8sClient is created'
        client1 instanceof K8sClient
        client1.config.server == 'http://k8s-server'

        when: 'second call within cache interval'
        def client2 = executor.getClient()
        then: 'returns the same cached instance'
        client2.is(client1)

        when: 'call after cache expiration'
        sleep(150)
        def client3 = executor.getClient()
        then: 'a new K8sClient instance is created'
        client3 instanceof K8sClient
        !client3.is(client1)
    }

}
