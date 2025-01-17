/*
 * Copyright 2021, Microsoft Corp
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

package nextflow.cloud.azure.config

import groovy.transform.CompileStatic
import nextflow.Global
import nextflow.Session

/**
 * Model Azure settings defined in the nextflow.config file
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzConfig {

    private AzCopyOpts azcopyOpts

    private AzStorageOpts storageOpts

    private AzBatchOpts batchOpts

    private AzRegistryOpts registryOpts

    private AzRetryConfig retryConfig

    private AzActiveDirectoryOpts activeDirectoryOpts

    private AzManagedIdentityOpts managedIdentityOpts

    AzConfig(Map azure) {
        this.batchOpts = new AzBatchOpts( (Map)azure.batch ?: Collections.emptyMap() )
        this.storageOpts = new AzStorageOpts( (Map)azure.storage ?: Collections.emptyMap() )
        this.registryOpts = new AzRegistryOpts( (Map)azure.registry ?: Collections.emptyMap() )
        this.azcopyOpts = new AzCopyOpts( (Map)azure.azcopy ?: Collections.emptyMap() )
        this.retryConfig = new AzRetryConfig( (Map)azure.retryPolicy ?: Collections.emptyMap() )
        this.activeDirectoryOpts = new AzActiveDirectoryOpts((Map) azure.activeDirectory ?: Collections.emptyMap())
        this.managedIdentityOpts = new AzManagedIdentityOpts((Map) azure.managedIdentity ?: Collections.emptyMap())
    }

    AzCopyOpts azcopy() { azcopyOpts }

    AzBatchOpts batch() { batchOpts }

    AzStorageOpts storage() { storageOpts }

    AzRegistryOpts registry() { registryOpts }

    AzRetryConfig retryConfig() { retryConfig }

    AzActiveDirectoryOpts activeDirectory() { activeDirectoryOpts }

    AzManagedIdentityOpts managedIdentity() { managedIdentityOpts }

    static AzConfig getConfig(Session session) {
        if( !session )
            throw new IllegalStateException("Missing Nextflow session")

        new AzConfig( (Map)session.config.azure ?: Collections.emptyMap()  )
    }

    static AzConfig getConfig() {
        getConfig(Global.session as Session)
    }
}
