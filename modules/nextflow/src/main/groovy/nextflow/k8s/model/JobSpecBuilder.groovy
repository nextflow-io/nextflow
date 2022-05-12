/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.k8s.model

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Object build for a K8s job specification
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Slf4j
class JobSpecBuilder {

    String jobName

    String namespace

    String workDir

    Integer backoffLimit = 0

    Map podSpec

    JobSpecBuilder withNamespace(String name) {
        this.namespace = name
        return  this
    }

    JobSpecBuilder withJobName(String name) {
        this.jobName = name
        return this
    }

    JobSpecBuilder withJobOptions(PodOptions opts) {
        if ( opts.backoffLimit != -1)
            backoffLimit = opts.backoffLimit
        return this
    }

    JobSpecBuilder withPodSpec(Map podSpec) {
        this.podSpec = podSpec
        return this
    }

    Map build() {
        assert this.jobName, 'Missing K8s jobName parameter'
        assert this.podSpec, 'Missing K8s podSpec parameter'

        final metadata = new LinkedHashMap<String,Object>()
        metadata.name = this.jobName
        metadata.namespace = this.namespace ?: 'default'

        /*final jobSpec = [
                template: podSpec,
                backoffLimit: backoffLimit
        ]*/

        final job = [
                apiVersion: 'batch/v1',
                kind: 'Job',
                metadata: metadata,
                spec: [ 
                   backoffLimit: backoffLimit,
                   template: podSpec
                ]
        ]
        return job
    }

}
