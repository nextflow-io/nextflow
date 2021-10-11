/*
 * Copyright 2021, Sage-Bionetworks
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

package nextflow.secret

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j

/**
 * Model a context to access secret values in nextflow config files
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class SecretsContext {

    SecretsContext() {}

    @Override
    Object getProperty(String name) {
        if( metaClass.hasProperty(name) )
            return metaClass.getProperty(this,name)
        else {
            return new SecretHolder(name)
        }
    }
}
