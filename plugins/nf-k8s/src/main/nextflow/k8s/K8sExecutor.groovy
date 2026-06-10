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

import com.google.common.cache.Cache
import com.google.common.cache.CacheBuilder
import groovy.transform.CompileStatic
import groovy.transform.Memoized
import groovy.util.logging.Slf4j
import nextflow.executor.Executor
import nextflow.fusion.FusionHelper
import nextflow.k8s.client.K8sClient
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.ServiceName
import org.pf4j.ExtensionPoint

/**
 * Implement the Kubernetes executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@ServiceName('k8s')
class K8sExecutor extends Executor implements ExtensionPoint {

    /**
     * Cache for the Kubernetes HTTP client. The client is refreshed periodically
     * so that the service account token is re-read when it expires.
     */
    private Cache<String, K8sClient> clientCache

    /**
     * @return The Kubernetes HTTP client. Delegates to a Guava cache that refreshes
     * the client (including the service account token) when the configured interval expires.
     */
    protected K8sClient getClient() {
        clientCache.get('client', () -> new K8sClient(k8sConfig.getClient()))
    }

    /**
     * @return The `k8s` configuration scope in the nextflow configuration object
     */
    @Memoized
    protected K8sConfig getK8sConfig() {
        new K8sConfig( (Map<String,Object>)session.config.k8s )
    }

    /**
     * Initialise the executor setting-up the kubernetes client configuration
     */
    @Override
    protected void register() {
        super.register()
        final k8sConfig = getK8sConfig()
        final refreshInterval = k8sConfig.clientRefreshInterval
        this.clientCache = CacheBuilder.newBuilder()
            .expireAfterWrite(refreshInterval.toMillis(), TimeUnit.MILLISECONDS)
            .build()
        final client = getClient()
        log.debug "[K8s] config=$k8sConfig; API client config=$client.config"
    }

    /**
     * @return {@code true} since containerised execution is managed by Kubernetes
     */
    boolean isContainerNative() {
        return true
    }

    @Override
    String containerConfigEngine() {
        return 'docker'
    }

    /**
     * @return A {@link TaskMonitor} associated to this executor type
     */
    @Override
    protected TaskMonitor createTaskMonitor() {
        TaskPollingMonitor.create(session, config, name, 100, Duration.of('5 sec'))
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

    @Override
    boolean isFusionEnabled() {
        return FusionHelper.isFusionEnabled(session)
    }
}
