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

import java.util.regex.Matcher
import java.util.regex.Pattern

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
    // TODO FIX THESE
    //static final private Pattern VERSION_JSON = Pattern.compile(".*/v(\\d+(?:\\.\\w+)*)-(\\w*)\\.json$")
    //static final private Pattern VERSION_TARGZ = Pattern.compile(".*/pkg\\/(\\d+(?:\\/\\w+)+)\\/fusion-(\\w+)\\.tar\\.gz$")
    final private Pattern VERSION_JSON = Pattern.compile('NEED-TO-FIX')
    final private Pattern VERSION_TARGZ = Pattern.compile('NEED-TO-FIX')

    boolean enabled
    String version

    FusionMetadata(Session session) {
        if( session.config.fusion as Map ) {
            final Map fusionConfig = session.config.fusion as Map
            this.enabled = fusionConfig.enabled as boolean
            // TODO work in progress
            this.version = this.enabled ? retrieveFusionVersion(fusionConfig) : null
        } else {
            this.enabled = false
            this.version = null
        }
    }

    FusionMetadata(Boolean enabled, String version) {
        this.enabled = enabled
        this.version = version
    }

    private String retrieveFusionVersion(Map config) {
        final String url = config.containerConfigUrl as String
        if( url && !url.isEmpty() && url.startsWith("https://fusionfs.seqera.io/") ) {
            final Matcher matcher_json = VERSION_JSON.matcher(url)
            if( matcher_json.matches() )
                return matcher_json.group(1)
            final Matcher matcher_targz = VERSION_TARGZ.matcher(url)
            if( matcher_targz.matches() )
                return matcher_targz.group(1).replaceAll("/", ".")
            return url
        }
        return null
    }

}