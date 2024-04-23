/*
 * Copyright 2013-2024, Seqera Labs
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
 *
 */

package nextflow.script

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.util.logging.Slf4j
import nextflow.Session
/**
 * Models Fusion metadata for Nextflow execution
 * 
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
@Slf4j
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode
class FusionMetadata {
    boolean enabled
    String version

    FusionMetadata(Session session) {
        if( session.config.fusion as Map ) {
            final Map fusionConfig = session.config.fusion as Map
            this.enabled = fusionConfig.enabled as boolean
        } else {
            this.enabled = false
        }
        // TODO work in progress
        this.version = 'not-defined-yet'
    }

}