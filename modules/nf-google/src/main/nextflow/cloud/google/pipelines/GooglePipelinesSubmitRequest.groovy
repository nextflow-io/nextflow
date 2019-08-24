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

package nextflow.cloud.google.pipelines

import com.google.api.services.genomics.v2alpha1.model.Mount
import groovy.transform.CompileStatic
import groovy.transform.ToString

/**
 * Models Google pipeline request for a Nextflow task executor
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
class GooglePipelinesSubmitRequest {

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

    String stagingScript

    String mainScript

    String unstagingScript

    Mount sharedMount

}
