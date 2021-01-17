/*
 * Copyright 2020, Microsoft Corp
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

import com.azure.storage.blob.nio.AzureFileSystem
import groovy.transform.CompileStatic
/**
 * Parse Azure settings from nextflow config file
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzStorageOpts {

    String accountKey
    String fileStores
    String accountName
    String sasToken

    AzStorageOpts(Map config) {
        assert config!=null
        this.accountKey = config.accountKey
        this.accountName = config.accountName
        this.fileStores = config.fileStores
        this.sasToken = config.sasToken
    }

    Map<String,Object> getEnv() {
        Map<String, Object> props = new HashMap<>();
        props.put(AzureFileSystem.AZURE_STORAGE_ACCOUNT_KEY, accountKey);
        props.put(AzureFileSystem.AZURE_STORAGE_FILE_STORES, fileStores);
        return props
    }

}
