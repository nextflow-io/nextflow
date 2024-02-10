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

package nextflow.container

import groovy.transform.CompileStatic

/**
 * Handle caching of remote Apptainer images
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class ApptainerCache extends SingularityCache {

    /** only for testing */
    protected ApptainerCache() {}

    ApptainerCache(ContainerConfig config, Map<String,String> env=null) {
        super(config, env)
    }

    @Override
    protected String getBinaryName() { 'apptainer' }
}
