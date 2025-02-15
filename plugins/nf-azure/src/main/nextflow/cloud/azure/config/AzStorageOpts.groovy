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
import nextflow.SysEnv
import nextflow.cloud.azure.batch.AzHelper
import nextflow.cloud.azure.nio.AzFileSystemProvider
import nextflow.util.Duration
/**
 * Parse Azure settings from nextflow config file
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzStorageOpts {

    String accountKey
    String accountName
    String sasToken
    Duration tokenDuration
    Map<String,AzFileShareOpts> fileShares


    AzStorageOpts(Map config, Map<String,String> env=SysEnv.get()) {
        assert config!=null
        this.accountKey = config.accountKey ?: env.get('AZURE_STORAGE_ACCOUNT_KEY')
        this.accountName = config.accountName ?: env.get('AZURE_STORAGE_ACCOUNT_NAME')
        this.sasToken = config.sasToken ?: env.get('AZURE_STORAGE_SAS_TOKEN')
        this.tokenDuration = (config.tokenDuration as Duration) ?: Duration.of('48h')
        this.fileShares = parseFileShares(config.fileShares instanceof Map ? config.fileShares as Map<String, Map>
                : Collections.<String,Map> emptyMap())

    }

    Map<String,Object> getEnv() {
        Map<String, Object> props = new HashMap<>();
        props.put(AzFileSystemProvider.AZURE_STORAGE_ACCOUNT_KEY, accountKey)
        props.put(AzFileSystemProvider.AZURE_STORAGE_ACCOUNT_NAME, accountName)
        props.put(AzFileSystemProvider.AZURE_STORAGE_SAS_TOKEN, sasToken)
        return props
    }

    static Map<String,AzFileShareOpts> parseFileShares(Map<String,Map> shares) {
        final result = new LinkedHashMap<String,AzFileShareOpts>()
        shares.each { Map.Entry<String, Map> entry ->
            result[entry.key] = new AzFileShareOpts(entry.value)
        }
        return result
    }
}
