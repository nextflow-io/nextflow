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
import nextflow.Global
import nextflow.Session
import nextflow.exception.AbortOperationException
import nextflow.plugin.Priority
/**
 * Implements the cloud cache factory
 *
 * @see CloudCacheStore
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
@Priority(-10)
class CloudCacheFactory extends CacheFactory {

    @Override
    protected CacheDB newInstance(UUID uniqueId, String runName, Path home) {
        if( !uniqueId ) throw new AbortOperationException("Missing cache `uuid`")
        if( !runName ) throw new AbortOperationException("Missing cache `runName`")
        final path = (Global.session as Session).cloudCachePath
        if( !path )
            throw new IllegalArgumentException("Cloud-cache path not defined - use either -with-cloudcache run option or NXF_CLOUDCACHE_PATH environment variable")
        final store = new CloudCacheStore(uniqueId, runName, path)
        return new CacheDB(store)
    }

}
