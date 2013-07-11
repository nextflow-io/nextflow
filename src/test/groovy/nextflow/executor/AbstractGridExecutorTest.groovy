/*
 * Copyright (c) 2012, the authors.
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

package nextflow.executor

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AbstractGridExecutorTest extends Specification {

    def 'test changeToScratchDir' () {

        setup:
        def executor = new AbstractGridExecutor() {
            @Override
            protected List<String> getSubmitCommandLine(TaskRun task) {
                return null  //To change body of implemented methods use File | Settings | File Templates.
            }
        }

        when:
        executor.taskConfig = new TaskConfig( )
        then:
        executor.changeToScratchDirectory() == null

        when:
        executor.taskConfig = new TaskConfig( [scratch: true] )
        then:
        executor.changeToScratchDirectory() == 'NF_SCRATCH=${TMPDIR:-`mktemp -d`} && cd $NF_SCRATCH'

        when:
        executor.taskConfig = new TaskConfig( [scratch: '$SOME_DIR'] )
        then:
        executor.changeToScratchDirectory() == 'NF_SCRATCH=$SOME_DIR && cd $NF_SCRATCH'

        when:
        executor.taskConfig = new TaskConfig( [scratch: '$SOME_DIR'] )
        then:
        executor.changeToScratchDirectory() == 'NF_SCRATCH=$SOME_DIR && cd $NF_SCRATCH'

        when:
        executor.taskConfig = new TaskConfig( [scratch: '/my/temp'] )
        then:
        executor.changeToScratchDirectory() == 'NF_SCRATCH=$(mktemp -d -p /my/temp) && cd $NF_SCRATCH'

        when:
        executor.taskConfig = new TaskConfig( [scratch: '/my/temp'] )
        then:
        executor.changeToScratchDirectory() == 'NF_SCRATCH=$(mktemp -d -p /my/temp) && cd $NF_SCRATCH'

    }

}
