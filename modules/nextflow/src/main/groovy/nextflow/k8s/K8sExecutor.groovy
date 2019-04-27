/*
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

package nextflow.k8s

import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.transform.PackageScope
import groovy.util.logging.Slf4j
import nextflow.executor.Executor
import nextflow.k8s.client.K8sClient
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.ServiceName
/**
 * Implement the Kubernetes executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@ServiceName('k8s')
class K8sExecutor extends Executor {

    /**
     * The Kubernetes HTTP client
     */
    private K8sClient client

    @PackageScope K8sClient getClient() {
        client
    }

    /**
     * @return The `k8s` configuration scope in the nextflow configuration object
     */
    @Memoized
    @PackageScope
    K8sConfig getK8sConfig() {
        new K8sConfig( (Map<String,Object>)session.config.k8s )
    }

    /**
     * Initialise the executor setting-up the kubernetes client configuration
     */
    @Override
    protected void register() {
        super.register()
        final config = k8sConfig.getClient()
        client = new K8sClient(config)
        log.debug "[K8s] API client config=$config"
    }

    /**
     * @return {@code true} since containerised execution is managed by Kubernetes
     */
    boolean isContainerNative() {
        return true
    }

    /**
     * @return A {@link TaskMonitor} associated to this executor type
     */
    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, name, 100, Duration.of('5 sec'))
    }

    /**
     * Creates a {@link TaskHandler} for the given {@link TaskRun} instance
     *
     * @param task A {@link TaskRun} instance representing a process task to be executed
     * @return A {@link K8sTaskHandler} instance modeling the execution in the K8s cluster
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir
        log.trace "[K8s] launching process > ${task.name} -- work folder: ${task.workDirStr}"
        new K8sTaskHandler(task,this)
    }
}
