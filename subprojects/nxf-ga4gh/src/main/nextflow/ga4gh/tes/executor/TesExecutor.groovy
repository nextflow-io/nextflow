/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.ga4gh.tes.executor

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.exception.AbortOperationException
import nextflow.executor.Executor
import nextflow.ga4gh.tes.client.ApiClient
import nextflow.ga4gh.tes.client.api.TaskServiceApi
import nextflow.processor.TaskHandler
import nextflow.processor.TaskMonitor
import nextflow.processor.TaskPollingMonitor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.ServiceName
/**
 * Experimental TES executor
 *
 * See https://github.com/ga4gh/task-execution-schemas/
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@ServiceName('tes')
class TesExecutor extends Executor {

    static private TaskServiceApi client

    @Override
    void register() {
        if( session.binDir && !session.binDir.empty() ) {
            session.abort()
            throw new AbortOperationException("ERROR: TES executor does not allow the use of custom scripts in the `bin` folder")
        }

        super.register()

        client = new TaskServiceApi( new ApiClient(basePath: getEndPoint()) )
    }

    protected String getDisplayName() {
        return "$name [${getEndPoint()}]"
    }

    TaskServiceApi getClient() {
        client
    }

    protected String getEndPoint() {
        def result = session.getConfigAttribute('executor.tes.endpoint', 'http://localhost:8000')
        log.debug "[TES] endpoint=$result"
        return result
    }

    /**
     * @return {@code true} whenever the containerization is managed by the executor itself
     */
    boolean isContainerNative() {
        return true
    }

    /**
     * Create a a queue holder for this executor
     *
     * @return
     */
    TaskMonitor createTaskMonitor() {
        return TaskPollingMonitor.create(session, name, 100, Duration.of('1 sec'))
    }


    /*
     * Prepare and launch the task in the underlying execution platform
     */
    @Override
    TaskHandler createTaskHandler(TaskRun task) {
        assert task
        assert task.workDir
        log.debug "[TES] Launching process > ${task.name} -- work folder: ${task.workDir}"
        new TesTaskHandler(task, this)
    }
}



