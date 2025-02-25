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

package nextflow.cache

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.exception.AbortOperationException
import nextflow.plugin.Priority

/**
 * Implements the default cache factory
 *
 * @see DefaultCacheStore
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@Priority(0)
class DefaultCacheFactory extends CacheFactory {

    @Override
    protected CacheDB newInstance(UUID uniqueId, String runName, Path home) {
        if( !uniqueId ) throw new AbortOperationException("Missing cache `uuid`")
        if( !runName ) throw new AbortOperationException("Missing cache `runName`")
        final store = new DefaultCacheStore(uniqueId, runName, home)
        return new CacheDB(store)
    }

}
