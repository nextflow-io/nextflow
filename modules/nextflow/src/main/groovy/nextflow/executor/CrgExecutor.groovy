/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.executor

import groovy.util.logging.Slf4j
import nextflow.processor.TaskRun
/**
 * An executor specialised for CRG cluster
 */
@Slf4j
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

        if( task.config.cpus>1 && !task.config.penv ) {
            log.debug 'Parallel environment not specified -- Using default value: `smp`'
            task.config.penv = 'smp'
        }

        super.getDirectives(task, result)

        if( task.config.getDisk() ) {
            result << "-l" << "disk=${task.config.getDisk().toMega()}M"
        }

        if( task.container && task.isDockerEnabled() ) {
            //  this will export the SGE_BINDING environment variable used to set Docker cpuset
            result << '-binding' << "env linear:${task.config.cpus}"

            // when it is a parallel job add 'reserve' flag
            if( task.config.cpus>1 ) {
                result << '-R' << 'y'
            }

            // request the docker image as a soft resource
            result << '-soft' << "-l docker_images=*;${task.container};*"
        }

        return result
    }


    @Override
    protected BashWrapperBuilder createBashWrapperBuilder(TaskRun task) {

        def builder = super.createBashWrapperBuilder(task)

        // When job is execute in a docker container
        // The Univa scheduler must allocate the required cores for the job execution
        // The variable '$SGE_BINDING' must contain the cores to be used
        if( task.container && task.isDockerEnabled() ) {
            def opt = task.containerConfig.legacy ? '--cpuset' : '--cpuset-cpus'
            def str = "\n"
            str += "cpuset=\${cpuset:=''}\n"
            str += "[[ \$SGE_BINDING ]] && cpuset=\"$opt \$(echo \$SGE_BINDING | sed 's/ /,/g')\"\n"
            builder.setContainerCpuset('$cpuset')
            builder.headerScript += str
        }

        return builder
    }

}
