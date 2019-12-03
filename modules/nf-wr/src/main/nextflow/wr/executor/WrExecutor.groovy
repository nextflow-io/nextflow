/*
 * Copyright 2019, Genome Research Limited
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

package nextflow.wr.executor

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import java.nio.file.Paths
import java.nio.file.Path
import com.upplication.s3fs.S3Path
import groovy.transform.PackageScope

import nextflow.wr.client.WrRestApi
import nextflow.wr.processor.WrMonitor

import nextflow.executor.Executor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.exception.AbortOperationException
import nextflow.util.ServiceName
import nextflow.extension.FilesEx

/**
 * Executor that schedules jobs using wr as a backend, avoiding storage of
 * state on disk.
 *
 * See https://github.com/VertebrateResequencing/wr
 *
 * @author Sendu Bala <sb10@sanger.ac.uk>
 * Based on TesExecutor by Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@ServiceName('wr')
class WrExecutor extends Executor {

    static private WrRestApi client
    static private String managerDir
    static private String deployment
    static private String endpoint
    private Path remoteBinDir = null

    @Override
    void register() {
        super.register()

        // upload local binaries
        def disableBinDir = session.getExecConfigProp(name, 'disableRemoteBinDir', false)
        if (workDir instanceof S3Path && session.binDir && !session.binDir.empty() && !disableBinDir ) {
            def s3 = getTempDir()
            log.info "Uploading local `bin` scripts folder to ${s3.toUriString()}/bin"
            remoteBinDir = FilesEx.copyTo(session.binDir, s3)
        }
    }

    WrRestApi getClient() {
        if (!client) {
            resolveDeployment()
            resolveManagerDir()
            endpoint = getEndPoint()
            client = new WrRestApi(endpoint, getToken(), getCacertPath())
        }
        client
    }

    protected String getDisplayName() {
        return "$name [$endpoint]"
    }

    protected void resolveDeployment() {
        deployment = session.getConfigAttribute('executor.wr.deployment', 'production')
        log.debug "[wr] deployment=$deployment"
    }

    protected String getHomeDir() {
        System.getProperty('user.home')
    }

    protected String getManagerDirBaseName() {
        ".wr_$deployment"
    }

    protected void resolveManagerDir() {
        managerDir = Paths.get(getHomeDir(), getManagerDirBaseName())
    }

    protected Integer getUID() {
        ["id", "-u"].execute().text.trim() as Integer // *** is there a groovy/java built-in for getting uid?
    }

    protected Integer getPort(int uid) {
        // default port that wr listens on is 1021 + (uid * 4) + n
        // where n is 1 for production and 3 for development
        int port = 1021 + (uid * 4) + 1
        if (deployment == 'development') {
            port += 2
        }
        return port
    }

    protected String getEndPoint() {
        int port = getPort(getUID())
        def result = session.getConfigAttribute('executor.wr.endpoint', "https://localhost:$port")
        log.debug "[wr] endpoint=$result"
        return result
    }

    protected String getToken() {
        String path = session.getConfigAttribute('executor.wr.tokenpath', Paths.get(managerDir, "client.token"))
        log.debug "[wr] tokenpath=$path"
        String result = new File(path).text
        return result
    }

    protected String getCacertPath() {
        String path = session.getConfigAttribute('executor.wr.cacertpath', Paths.get(managerDir, "ca.pem"))
        log.debug "[wr] cacertpath=$path"
        return path
    }

    @PackageScope
    Path getRemoteBinDir() {
        remoteBinDir
    }

    /**
     * @return {@code false} whenever the containerization is not managed by the executor itself
     */
    boolean isContainerNative() {
        return false
    }

    /**
     * Create a queue holder for this executor
     *
     * @return
     */
    TaskMonitor createTaskMonitor() {
        return WrMonitor.create(session, getClient())
    }

    /*
     * Prepare and launch the task in the underlying execution platform
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir
        log.debug "[wr] Launching process > ${task.name} -- work folder: ${task.workDir}"
        new WrTaskHandler(task, this)
    }
    
}
