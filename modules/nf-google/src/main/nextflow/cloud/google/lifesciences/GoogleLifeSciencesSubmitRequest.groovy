/*
 * Copyright 2019, Google Inc
 * Copyright 2018, WuxiNextcode
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

package nextflow.cloud.google.lifesciences

import java.nio.file.Path

import com.google.api.services.lifesciences.v2beta.model.Mount
import groovy.transform.CompileStatic
import groovy.transform.ToString
import nextflow.executor.res.AcceleratorResource
import nextflow.processor.TaskRun

/**
 * Models Google pipeline request for a Nextflow task executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
class GoogleLifeSciencesSubmitRequest implements GoogleLifeSciencesTaskDirWrangler {

    String machineType

    String project

    List<String> zone

    List<String> region

    String diskName

    Integer diskSizeGb

    boolean preemptible

    String taskName

    String containerImage

    String fileCopyImage

    Mount sharedMount

    AcceleratorResource accelerator

    String location

    Path workDir


    String getStagingScript() {
        String result = "set -x; "
        result += '{ '
        result += "cd ${localTaskDir}; "
        result += "gsutil -m -q cp $remoteTaskDir/${TaskRun.CMD_RUN} .; "
        result += "bash ${TaskRun.CMD_RUN} nxf_stage; "
        result += "} 2>&1 > $localTaskDir/${TaskRun.CMD_LOG}"
        return result
    }

    String getMainScript() {
        "{ cd $localTaskDir; bash ${TaskRun.CMD_RUN}; } 2>&1 | tee -a $localTaskDir/${TaskRun.CMD_LOG}"
    }

    String getUnstagingScript() {
        def result = "set -x; "
        result += "{ cd $localTaskDir; bash ${TaskRun.CMD_RUN} nxf_unstage; } 2>&1 | tee -a $localTaskDir/${TaskRun.CMD_LOG}; "
        result += "gsutil -m -q cp -R $localTaskDir/${TaskRun.CMD_LOG} ${remoteTaskDir}/${TaskRun.CMD_LOG} || true; "
        result += "[[ \$GOOGLE_LAST_EXIT_STATUS -gt 0 || \$NXF_DEBUG -gt 0 ]] && { gsutil -m -q cp -R /google/ ${remoteTaskDir}; } || rm -rf $localTaskDir"
        return result
    }
}
