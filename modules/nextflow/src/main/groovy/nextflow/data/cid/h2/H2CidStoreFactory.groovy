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

package nextflow.data.cid.h2

import groovy.transform.CompileStatic
import nextflow.data.cid.CidStore
import nextflow.data.cid.CidStoreFactory
import nextflow.data.config.DataConfig

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class H2CidStoreFactory extends CidStoreFactory {

    @Override
    boolean canOpen(DataConfig config) {
        return config.store.location.startsWith('jdbc:h2:')
    }

    @Override
    protected CidStore newInstance(DataConfig config) {
        return new H2CidStore()
    }
}
