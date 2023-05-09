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

package nextflow.k8s.client

import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.k8s.client.res.JobProvider
import nextflow.k8s.client.res.MPIJobProvider
import nextflow.k8s.client.res.PodProvider
/**
 * Kubernetes API client
 *
 * Tip: use the following command to find out your kubernetes master node
 *   kubectl cluster-info
 *
 * See
 *   https://v1-8.docs.kubernetes.io/docs/api-reference/v1.8/#-strong-api-overview-strong-
 *
 * Useful cheatsheet
 *   https://kubernetes.io/docs/reference/kubectl/cheatsheet/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class K8sClient {

    @Delegate
    private K8sRequestDelegate delegate

    private Map<String,ResourceProvider> resourceProviders

    K8sClient() {
        this(new ClientConfig())
    }

    ClientConfig getConfig() { delegate.config }

    /**
     * Creates a Kubernetes client using the configuration setting provided by the specified
     * {@link ClientConfig} instance.
     *
     * @param config
     */
    K8sClient(ClientConfig config) {
        this.delegate = new K8sRequestDelegate(config)
        this.resourceProviders = [:]
        registerProvider('Pod', PodProvider.class)
        registerProvider('Job', JobProvider.class)
        registerProvider('MPIJob', MPIJobProvider.class)
    }

    K8sClient(ClientConfig config, Map resourceProviders) {
        this.delegate = new K8sRequestDelegate(config)
        this.resourceProviders = resourceProviders
    }

    void registerProvider(String name, Class<? extends ResourceProvider> clazz) {
        final provider = clazz.newInstance().withDelegate(delegate)
        this.resourceProviders.put(name, provider)
    }

    K8sResponseJson create(String type, Map spec, Path saveYamlPath=null) {
        resourceProviders[type].create(spec, saveYamlPath)
    }

    K8sResponseJson delete(String type, String name) {
        resourceProviders[type].delete(name)
    }

    K8sResponseJson list(String type, boolean allNamespaces=false) {
        resourceProviders[type].list(allNamespaces)
    }

    InputStream fetchLog(Map params=[:], String type, String name) {
        resourceProviders[type].fetchLog(params, name)
    }

    Map getState( String type, String name ) {
        resourceProviders[type].getState(name)
    }

    String podHostname(String name){
        ((PodProvider)resourceProviders['Pod']).getHostname(name)
    }

    K8sResponseJson configCreate(String name, Map data) {

        final spec = [
                apiVersion: 'v1',
                kind: 'ConfigMap',
                metadata: [ name: name, namespace: config.namespace ],
                data: data
        ]

        final action = "/api/v1/namespaces/${config.namespace}/configmaps"
        final body = JsonOutput.toJson(spec)
        def resp = makeRequest('POST', action, body)
        trace('POST', action, resp.text)
        return new K8sResponseJson(resp.text)
    }

    K8sResponseJson configDelete(String name) {
        final action = "/api/v1/namespaces/${config.namespace}/configmaps/${name}"
        def resp = makeRequest('DELETE', action)
        trace('DELETE', action, resp.text)
        return new K8sResponseJson(resp.text)
    }

    K8sResponseJson configDeleteAll() {
        final action = "/api/v1/namespaces/${config.namespace}/configmaps"
        def resp = makeRequest('DELETE', action)
        trace('DELETE', action, resp.text)
        return new K8sResponseJson(resp.text)
    }

    K8sResponseJson secretDescribe(String name) {
        assert name
        final action = "/api/v1/namespaces/${config.namespace}/secrets/${name}"
        final resp = makeRequest('GET', action)
        trace('GET', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    K8sResponseJson secretList() {
        final action = "/api/v1/namespaces/${config.namespace}/secrets"
        final resp = makeRequest('GET', action)
        trace('GET', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    K8sResponseJson volumeClaimRead(String name) {
        final action = "/api/v1/namespaces/${config.namespace}/persistentvolumeclaims/${name}"
        def resp = makeRequest('GET', action)
        trace('GET', action, resp.text)
        return new K8sResponseJson(resp.text)
    }


}

