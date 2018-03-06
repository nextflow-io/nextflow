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

package nextflow.ga4gh.tes.server.verticle

import nextflow.Session
import nextflow.executor.Executor
import nextflow.processor.ProcessConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.script.BaseScript
import nextflow.script.TaskBody
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TesTaskProcessorAdaptor extends TaskProcessor {

    TesTaskProcessorAdaptor(Session session, Executor executor) {
        this.session = session
        this.executor = executor
        this.ownerScript = new BaseScript() { @Override Object run() { return null } }
        this.config = new ProcessConfig(ownerScript)
    }


    @Override
    TaskBody getTaskBody() {
        return new TaskBody(null,null, ) {
            @Override Set<String> getValNames() { Collections.emptySet() }
        }
    }

    @Override
    protected finalizeTask(TaskRun task) {
        // disable task finalisation
    }
}
