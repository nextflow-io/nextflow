/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

import groovy.transform.PackageScope
import nextflow.processor.TaskRun
/**
 * An executor specialised for CRG cluster
 */
class CrgExecutor extends SgeExecutor {

    /**
     * The jobs directives used for the job submission
     *
     * @param task The {@link TaskRun} instance to be submitted
     * @param result A {@link List} object to which are added the directive tokens
     * @return The a list of string containing the directives names and values
     */
    @Override
    List<String> getDirectives(TaskRun task, List<String> result) {

        if( taskConfig.cpus>1 && !taskConfig.penv ) {
            log.debug 'Parallel environment not specified -- Using default value: `smp`'
            taskConfig.penv = 'smp'
        }

        super.getDirectives(task, result)

        if( task.container && isDockerEnabled() ) {
            result << '-soft' << '-l' << "docker_images=${task.container}"
        }

        return result
    }

    /**
     * @return The value of the {@code docker.enabled} configuration setting defined in the
     *  nextflow.config file
     */
    @PackageScope
    boolean isDockerEnabled() {
        Map dockerConf = session.config.docker as Map
        dockerConf?.enabled?.toString() == 'true'
    }


}
