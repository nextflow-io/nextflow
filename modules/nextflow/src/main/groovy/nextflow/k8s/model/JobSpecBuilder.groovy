/*
 * Copyright 2013-2023, Seqera Labs
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
 * Builder for a K8s Job specification
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Slf4j
class JobSpecBuilder extends PodSpecBuilder {

    @Override
    Map build() {
        final metadata = super.buildMetadata()
        final spec = [
            backoffLimit: 0,
            template: [
                spec: super.buildSpec()
            ]
        ]

        return [
            apiVersion: 'batch/v1',
            kind: 'Job',
            metadata: metadata,
            spec: spec
        ]
    }

}
