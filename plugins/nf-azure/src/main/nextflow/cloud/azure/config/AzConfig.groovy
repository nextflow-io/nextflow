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
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName
import nextflow.script.dsl.Description

/**
 * Model Azure settings defined in the nextflow.config file
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ScopeName("azure")
@Description("""
    The `azure` scope allows you to configure the interactions with Azure, including Azure Batch and Azure Blob Storage.
""")
@CompileStatic
class AzConfig implements ConfigScope {

    private AzCopyOpts azcopy

    private AzStorageOpts storage

    private AzBatchOpts batch

    private AzRegistryOpts registry

    private AzRetryConfig retryPolicy

    private AzActiveDirectoryOpts activeDirectory

    private AzManagedIdentityOpts managedIdentity

    /* required by extension point -- do not remove */
    AzConfig() {}

    AzConfig(Map azure) {
        this.batch = new AzBatchOpts( (Map)azure.batch ?: Collections.emptyMap() )
        this.storage = new AzStorageOpts( (Map)azure.storage ?: Collections.emptyMap() )
        this.registry = new AzRegistryOpts( (Map)azure.registry ?: Collections.emptyMap() )
        this.azcopy = new AzCopyOpts( (Map)azure.azcopy ?: Collections.emptyMap() )
        this.retryPolicy = new AzRetryConfig( (Map)azure.retryPolicy ?: Collections.emptyMap() )
        this.activeDirectory = new AzActiveDirectoryOpts((Map) azure.activeDirectory ?: Collections.emptyMap())
        this.managedIdentity = new AzManagedIdentityOpts((Map) azure.managedIdentity ?: Collections.emptyMap())
    }

    AzCopyOpts azcopy() { azcopy }

    AzBatchOpts batch() { batch }

    AzStorageOpts storage() { storage }

    AzRegistryOpts registry() { registry }

    AzRetryConfig retryConfig() { retryPolicy }

    AzActiveDirectoryOpts activeDirectory() { activeDirectory }

    AzManagedIdentityOpts managedIdentity() { managedIdentity }

    static AzConfig getConfig(Session session) {
        if( !session )
            throw new IllegalStateException("Missing Nextflow session")

        new AzConfig( (Map)session.config.azure ?: Collections.emptyMap()  )
    }

    static AzConfig getConfig() {
        getConfig(Global.session as Session)
    }
}
