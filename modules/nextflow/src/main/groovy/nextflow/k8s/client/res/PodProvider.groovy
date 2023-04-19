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
import nextflow.exception.K8sOutOfCpuException
import nextflow.exception.K8sOutOfMemoryException
import nextflow.exception.NodeTerminationException
import nextflow.exception.ProcessFailedException
import nextflow.k8s.client.ResourceProvider
import nextflow.k8s.client.K8sRequestDelegate
import nextflow.k8s.client.K8sResponseException
import nextflow.k8s.client.K8sResponseJson
import nextflow.k8s.client.PodUnschedulableException
import org.yaml.snakeyaml.Yaml
/**
 * Models API operations for Pods.
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class PodProvider implements ResourceProvider {

    @Delegate
    protected K8sRequestDelegate delegate

    @Override
    ResourceProvider withDelegate(K8sRequestDelegate delegate) {
        this.delegate = delegate
        return this
    }

    /**
     * Create a pod.
     *
     * See:
     *   https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#create-pod-v1-core
     *
     * @param spec
     * @param saveYamlPath
     */
    @Override
    K8sResponseJson create(Map spec, Path saveYamlPath=null) {

        log.debug "PodProvider::create(spec=$spec, saveYamlPath=$saveYamlPath)"

        if( saveYamlPath ) try {
            saveYamlPath.text = new Yaml().dump(spec).toString()
        }
        catch( Exception e ) {
            log.debug "WARN: unable to save request yaml -- cause: ${e.message ?: e}"
        }

        final req = JsonOutput.toJson(spec)
        final action = "/api/v1/namespaces/${config.namespace}/pods"
        final resp = makeRequest('POST', action, req)
        trace('POST', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    /**
     * Delete a pod.
     *
     * See:
     *   https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#delete-pod-v1-core
     *
     * @param name
     */
    @Override
    K8sResponseJson delete(String name) {
        final action = "/api/v1/namespaces/${config.namespace}/pods/${name}"
        final resp = makeRequest('DELETE', action)
        trace('DELETE', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    /**
     * List all pods in a namespace (or all namespaces).
     *
     * See:
     *   https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#list-pod-v1-core
     *   https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#list-all-namespaces-pod-v1-core
     *
     * @param allNamespaces
     */
    @Override
    K8sResponseJson list(boolean allNamespaces=false) {
        final String action = allNamespaces ? "pods" : "namespaces/${config.namespace}/pods"
        final resp = makeRequest('GET', "/api/v1/${action}")
        trace('GET', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    /**
     * Fetch the log of a pod.
     *
     * See:
     *   https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#read-log-pod-v1-core
     *
     * @param params
     * @param name
     */
    @Override
    InputStream fetchLog(Map params=[:], String name) {
        // -- compose the request action uri
        def action = "/api/v1/namespaces/${config.namespace}/pods/${name}/log"
        def count = 0
        for( String key : params.keySet() )
            action += "${count++ == 0 ? '?' : '&'}${key}=${params.get(key)}"
        // -- submit request
        def resp = makeRequest('GET', action)
        resp.stream
    }

    /**
     * Get the state of a pod.
     *
     * @param name
     * @return
     *     A {@link Map} representing the container state object as shown below:
     *     <code>
     *     {
     *         "terminated": {
     *             "exitCode": 127,
     *             "reason": "ContainerCannotRun",
     *             "message": "OCI runtime create failed: container_linux.go:296: starting container process caused \"exec: \\\"bash\\\": executable file not found in $PATH\": unknown",
     *             "startedAt": "2018-01-12T22:04:25Z",
     *             "finishedAt": "2018-01-12T22:04:25Z",
     *             "containerID": "docker://730ef2e05be72ffc354f2682b4e8300610812137b9037b726c21e5c4e41b6dda"
     *         }
     *     }
     *     </code>
     *
     *     An empty map is returned if the pod status is `Pending` and the container status is not
     *     yet available.
     *
     *     See:
     *       https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#read-status-pod-v1-core
     */
    @Override
    Map getState(String name) {
        final resp = getPodStatus(name)
        final status = resp.status as Map
        final containerStatuses = status?.containerStatuses as List<Map>

        if( containerStatuses?.size() > 0 ) {
            final container = containerStatuses.get(0)
            // note: when the pod is created by a Job submission
            // the `name` does not match the container name because it
            // contains a suffix random generated by K8s pod scheduler
            if( !container.name || !name.startsWith(container.name.toString()) )
                throw new K8sResponseException("K8s invalid status for pod: ${name} (unexpected container name: ${container.name})", resp)

            if( !container.state )
                throw new K8sResponseException("K8s invalid status for pod: ${name} (missing state object)", resp)

            final state = container.state as Map
            if( state.waiting instanceof Map ) {
                def waiting = state.waiting as Map
                checkInvalidWaitingState(waiting, resp)
            }
            return state
        }

        if( status?.phase == 'Pending' ) {
            if( status.conditions instanceof List ) {
                final allConditions = status.conditions as List<Map>
                final cond = allConditions.find { cond -> cond.type == 'PodScheduled' }
                if( cond.reason == 'Unschedulable' ) {
                    def message = "K8s pod cannot be scheduled"
                    if( cond.message ) message += " -- ${cond.message}"
                    //def cause = new K8sResponseException(resp)
                    log.warn1(message)
                }
            }
            // undetermined status -- return an empty response
            return Collections.emptyMap()
        }

        if( status?.phase == 'Failed' ) {
            def msg = "K8s pod '$name' execution failed"
            if( status.reason ) msg += " - reason: ${status.reason}"
            if( status.message ) msg += " - message: ${status.message}"
            switch ( status.reason ) {
                case 'OutOfcpu':    throw new K8sOutOfCpuException(msg)
                case 'OutOfmemory': throw new K8sOutOfMemoryException(msg)
                case 'Shutdown':    throw new NodeTerminationException(msg)
                default:            throw new ProcessFailedException(msg)
            }
        }

        throw new K8sResponseException("K8s undetermined status conditions for pod $name", resp)
    }

    protected K8sResponseJson getPodStatus(String name) {
        try {
            return getPodStatus0(name)
        }
        catch (K8sResponseException err) {
            if( err.response.code == 404 && isKindPods(err.response) ) {
                // this may happen when K8s node is shutdown and the pod is evicted
                // therefore process exception is thrown so that the failure
                // can be managed by the nextflow as retry-able execution
                throw new NodeTerminationException("Unable to find pod $name - The pod may be evicted by a node shutdown event")
            }
            throw err
        }
    }

    protected K8sResponseJson getPodStatus0(String name) {
        final action = "/api/v1/namespaces/${config.namespace}/pods/${name}/status"
        final resp = makeRequest('GET', action)
        trace('GET', action, resp.text)
        return new K8sResponseJson(resp.text)
    }

    protected boolean isKindPods(K8sResponseJson resp) {
        if( resp.details instanceof Map ) {
            final details = (Map) resp.details
            return details.kind == 'pods'
        }
        return false
    }

    protected void checkInvalidWaitingState( Map waiting, K8sResponseJson resp ) {
        if( waiting.reason == 'ErrImagePull' || waiting.reason == 'ImagePullBackOff') {
            def message = "K8s pod image cannot be pulled"
            if( waiting.message ) message += " -- ${waiting.message}"
            final cause = new K8sResponseException(resp)
            throw new PodUnschedulableException(message, cause)
        }
        if( waiting.reason == 'CreateContainerConfigError' ) {
            def message = "K8s pod configuration failed"
            if( waiting.message ) message += " -- ${waiting.message}"
            final cause = new K8sResponseException(resp)
            throw new PodUnschedulableException(message, cause)
        }
        if( waiting.reason =~ /.+Error$/ ) {
            def message = "K8s pod waiting on unknown error state"
            if( waiting.message ) message += " -- ${waiting.message}"
            final cause = new K8sResponseException(resp)
            throw new PodUnschedulableException(message, cause)
        }
        final status = resp.status as Map
        if( status?.phase == 'Failed' ) {
            def message = "K8s pod in Failed state"
            final cause = new K8sResponseException(resp)
            throw new PodUnschedulableException(message, cause)
        }
    }

    /**
     * Get the node hostname of a pod.
     *
     * @param name
     */
    String getHostname(String name){
        final K8sResponseJson resp = getPodStatus(name)
        (resp?.spec as Map)?.nodeName as String
    }

}
