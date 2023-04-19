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

package nextflow.k8s.client.res

import java.nio.file.Path

import groovy.json.JsonOutput
import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.k8s.client.K8sClient
import nextflow.k8s.client.K8sResponseException
import nextflow.k8s.client.K8sResponseJson
import org.yaml.snakeyaml.Yaml
/**
 * Models API operations for MPIJobs.
 *
 * See:
 *   https://github.com/kubeflow/mpi-operator
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class MPIJobProvider extends JobProvider {

    /**
     * Create an MPIJob.
     *
     * @param spec
     * @param saveYamlPath
     */
    @Override
    K8sResponseJson create(Map spec, Path saveYamlPath=null) {

        if( saveYamlPath ) try {
            saveYamlPath.text = new Yaml().dump(spec).toString()
        }
        catch( Exception e ) {
            log.debug "WARN: unable to save request yaml -- cause: ${e.message ?: e}"
        }

        final req = JsonOutput.toJson(spec)
        final action = "/apis/kubeflow.org/v2beta1/namespaces/${config.namespace}/mpijobs"
        final resp = makeRequest('POST', action, req)
        trace('POST', action, resp.text)
        return new K8sResponseJson(resp.text)
    }

    /**
     * Delete an MPIJob.
     *
     * @param name
     */
    @Override
    K8sResponseJson delete(String name) {
        final action1 = "/apis/kubeflow.org/v2beta1/namespaces/${config.namespace}/mpijobs/${name}"
        final resp1 = makeRequest('DELETE', action1)
        trace('DELETE', action1, resp1.text)
        new K8sResponseJson(resp1.text)
    }

    /**
     * List all MPIJobs in a namespace (or all namespaces).
     *
     * @param allNamespaces
     */
    @Override
    K8sResponseJson list(boolean allNamespaces=false) {
        final String action = allNamespaces ? "mpijobs" : "namespaces/${config.namespace}/mpijobs"
        final resp = makeRequest('GET', "/apis/kubeflow.org/v2beta1/${action}")
        trace('GET', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    /**
     * Get the state of an MPIJob.
     *
     * @param name
     */
    @Override
    Map getState(String name) {
        // get launcher pod name
        name = "${name}-launcher"

        // try to get the pod state multiple times with a backoff,
        // since the launcher pod may not be created yet
        for( int i = 1; i <= 5; i++ ) {
            try {
                return super.getState(name) 
            } 
            catch( K8sResponseException err ) {
                if( err.response.code != 404 )
                    throw err
                else {
                    final long delay = (Math.pow(3, i - 1) as long) * 250
                    sleep( delay )
                }
            }
        }
    }

    @Override
    InputStream fetchLog(Map params=[:], String name) {
        throw new UnsupportedOperationException()
    }

}

