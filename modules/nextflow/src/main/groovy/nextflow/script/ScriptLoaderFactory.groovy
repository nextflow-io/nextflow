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
import groovy.util.logging.Slf4j
import nextflow.Session
import nextflow.SysEnv
import nextflow.script.parser.v1.ScriptLoaderV1
import nextflow.script.parser.v2.ScriptLoaderV2

/**
 * Factory for creating an instance of {@link ScriptLoader}.
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@Slf4j
@CompileStatic
class ScriptLoaderFactory {

    static ScriptLoader create(Session session) {
        final parser = SysEnv.get('NXF_SYNTAX_PARSER', 'v1')
        if( parser == 'v1' ) {
            return new ScriptLoaderV1(session)
        }
        if( parser == 'v2' ) {
            log.debug "Using script parser v2"
            return new ScriptLoaderV2(session)
        }
        throw new IllegalStateException("Invalid NXF_SYNTAX_PARSER setting -- should be either 'v1' or 'v2'")
    }

}
