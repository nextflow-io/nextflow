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

import nextflow.wr.client.WrRestApi
import nextflow.wr.processor.WrMonitor

import nextflow.executor.Executor
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskRun
import nextflow.exception.AbortOperationException
import nextflow.util.ServiceName

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

    static private String token
    static private String defaultManagerDir
    static private String endpoint
    static private String cacertPath
    static private WrRestApi client

    @Override
    void register() {
        if( session.binDir && !session.binDir.empty() ) {
            session.abort()
            throw new AbortOperationException("ERROR: wr executor does not allow the use of custom scripts in the `bin` folder")
        }

        super.register()
    }

    protected String getDisplayName() {
        return "$name [$endpoint]"
    }

    WrRestApi getClient() {
        if (!client) {
            defaultManagerDir = Paths.get(System.getProperty('user.home'), ".wr_production")
            endpoint = getEndPoint()
            client = new WrRestApi(endpoint, getToken(), getCacertPath())
        }
        client
    }

    protected String getEndPoint() {
        // default port that wr listens on is 1021 + (uid * 4) + 1
        // *** note that this will probably only work on linux/mac os, but wr
        // probably only works fully on those as well...
        int uid = ["id", "-u"].execute().text.trim() as Integer
        int port = 1021 + (uid * 4) + 1

        def result = session.getConfigAttribute('executor.endpoint', "https://localhost:$port")
        log.debug "[wr] endpoint=$result"
        return result
    }

    protected String getToken() {
        String path = session.getConfigAttribute('executor.tokenpath', Paths.get(defaultManagerDir, "client.token"))
        log.debug "[wr] tokenpath=$path"
        String result = new File(path).text
        return result
    }

    protected String getCacertPath() {
        String path = session.getConfigAttribute('executor.cacertpath', Paths.get(defaultManagerDir, "ca.pem"))
        log.debug "[wr] cacertpath=$path"
        return path
    }

    /**
     * @return {@code false} whenever the containerization is managed by the executor itself
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
