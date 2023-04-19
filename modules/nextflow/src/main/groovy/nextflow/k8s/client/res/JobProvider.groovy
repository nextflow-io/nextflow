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
import nextflow.exception.NodeTerminationException
import nextflow.exception.ProcessFailedException
import nextflow.k8s.client.K8sResponseException
import nextflow.k8s.client.K8sResponseJson
import org.yaml.snakeyaml.Yaml
/**
 * Models API operations for Jobs.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class JobProvider extends PodProvider {

    /**
     * Create a job.
     *
     * See:
     *   https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#create-job-v1-batch
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
        final action = "/apis/batch/v1/namespaces/${config.namespace}/jobs"
        final resp = makeRequest('POST', action, req)
        trace('POST', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    /**
     * Delete a job.
     *
     * See:
     *   https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#delete-job-v1-batch
     *
     * @param name
     */
    @Override
    K8sResponseJson delete(String name) {
        // delete all pods managed by the job
        final action = "/api/v1/namespaces/${config.namespace}/pods?labelSelector=job-name=${name}"
        final resp = makeRequest('GET', action)
        trace('GET', action, resp.text)
        final podList = new K8sResponseJson(resp.text)

        if( podList.kind == "PodList" ) { 
            for( Map item : podList.items as List<Map> ) {
                try {
                    final metadata = item.metadata as Map
                    final podName = metadata.name as String
                    super.delete(podName)
                }
                catch( K8sResponseException err ) {
                    if( err.response.code == 404 )
                        log.debug("Unable to delete pod for job ${name}, pod already deleted")
                    else
                        throw err
                }
            }
        }

        // delete job
        final action1 = "/apis/batch/v1/namespaces/${config.namespace}/jobs/${name}"
        final resp1 = makeRequest('DELETE', action1)
        trace('DELETE', action1, resp1.text)
        new K8sResponseJson(resp1.text)
    }

    /**
     * List all jobs in a namespace (or all namespaces).
     *
     * See:
     *   https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#list-job-v1-batch
     *   https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#list-all-namespaces-job-v1-batch
     *
     * @param allNamespaces
     */
    @Override
    K8sResponseJson list(boolean allNamespaces=false) {
        final String action = allNamespaces ? "jobs" : "namespaces/${config.namespace}/jobs"
        final resp = makeRequest('GET', "/apis/batch/v1/${action}")
        trace('GET', action, resp.text)
        new K8sResponseJson(resp.text)
    }

    /**
     * Fetch the log of the pod most recently created by a job.
     *
     * See:
     *   https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#read-log-pod-v1-core
     *
     * @param params
     * @param name
     */
    @Override
    InputStream fetchLog(Map params=[:], String name) {
        super.fetchLog(params, getLatestPodName(name))
    }

    /**
     * Get the state of a job.
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
     *       https://kubernetes.io/docs/reference/generated/kubernetes-api/v1.27/#read-status-job-v1-batch
     */
    @Override
    Map getState(String name) {
        final podName = getLatestPodName(name)
        if( podName ) {
            try {
                return super.getState(podName)
            } 
            catch( NodeTerminationException err ) {
                // pod might be deleted by control plane just after getLatestPodName() call
                // so try fallback to jobState
                log.warn1 "Pod not found for Job ${name}, it may have been deleted by control plane"
                return getJobStatus(name)           
            }
        }
        else {
            return getJobStatus(name)
        }
    }

    /**
     * Get the name of the pod most recently created by a job.
     *
     * @param name
     */
    protected String getLatestPodName(String name) {
        final action = "/api/v1/namespaces/${config.namespace}/pods?labelSelector=job-name=${name}"
        final resp = makeRequest('GET', action)
        trace('GET', action, resp.text)
        final podList = new K8sResponseJson(resp.text)

        if( podList.kind != 'PodList' )
            return null

        def latestPodName = null
        def latestTimestamp = '0000-00-00T00:00:00Z'

        for( Map item : podList.items as List<Map> ) {
            final metadata = item.metadata as Map
            final timestamp = metadata.creationTimestamp as String

            if( timestamp > latestTimestamp ) {
                latestPodName = metadata.name
                latestTimestamp = timestamp
            }
        }

        return latestPodName
    }

    protected Map getJobStatus(String name) {
        final resp = getJobStatus0(name)
        final status = resp.status as Map
        if( status?.succeeded == 1 && status.conditions instanceof List ) {
            final allConditions = status.conditions as List<Map>
            final cond = allConditions.find { cond -> cond.type == 'Complete' }

            if( cond?.status == 'True' ) {
                log.warn1("Job ${name} already completed and Pod is deleted")
                return [
                    terminated: [
                        exitcode: 0,
                        reason: "Completed",
                        startedAt: status.startTime,
                        finishedAt: status.completionTime,
                    ]
                ]
            } else {
                throw new ProcessFailedException("K8s Job ${name} succeeded but does not have `Complete` status: ${allConditions}")
            }
        }

        if( status?.failed && ((int)status.failed) > 0 ) {
            String message = 'unknown'
            if( status.conditions instanceof List ) {
                final allConditions = status.conditions as List<Map>
                final cond = allConditions.find { cond -> cond.type == 'Failed' }
                message = cond?.message
            }
            throw new ProcessFailedException("K8s Job ${name} execution failed: ${message}")
        }

        log.warn1("Pod not found for Job ${name} -- it may not have been scheduled yet")
        return Collections.emptyMap()
    }

    protected K8sResponseJson getJobStatus0(String name) {
        final action = "/apis/batch/v1/namespaces/${config.namespace}/jobs/${name}/status"
        final resp = makeRequest('GET', action)
        trace('GET', action, resp.text)
        new K8sResponseJson(resp.text)
    }

}
