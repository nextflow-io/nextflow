/*
 * Copyright 2024-2025, Seqera Labs
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

package nextflow.config.spec

import nextflow.script.types.Duration
import nextflow.script.types.MemoryUnit
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class SpecNodeTest extends Specification {

    def 'should infer process config options from process directives' () {
        given:
        def scope = SpecNode.ROOT.children.get('process')

        expect:
        scope.children.get('clusterOptions').types as Set == [ String, List, String[] ] as Set
        scope.children.get('cpus').types == [ Integer ]
        scope.children.get('errorStrategy').types == [ String ]
        scope.children.get('executor').types == [ String ]
        scope.children.get('ext').types == [ Map ]
        scope.children.get('memory').types == [ MemoryUnit ]
        scope.children.get('publishDir').types as Set == [ Map, String, List ] as Set
        scope.children.get('resourceLimits').types == [ Map ]
        scope.children.get('time').types == [ Duration ]
    }

}
