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
 *
 */

package nextflow.cloud.aws.cache

import java.nio.file.Path

import groovy.transform.CompileStatic
import nextflow.cache.CacheDB
import nextflow.cache.CacheFactory
import nextflow.exception.AbortOperationException
/**
 * Implements the S3 cache factory
 *
 * @see S3CacheStore
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
@CompileStatic
class S3CacheFactory extends CacheFactory {

    @Override
    protected String getName() { 's3' }

    @Override
    protected CacheDB newInstance(UUID uniqueId, String runName, Path home) {
        if( !uniqueId ) throw new AbortOperationException("Missing cache `uuid`")
        if( !runName ) throw new AbortOperationException("Missing cache `runName`")
        final store = new S3CacheStore(uniqueId, runName, home)
        return new CacheDB(store)
    }

}
