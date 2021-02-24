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
import nextflow.cloud.azure.nio.AzFileSystemProvider
import nextflow.util.Duration

/**
 * Parse Azure settings from nextflow config file
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzStorageOpts {

    private Map<String,String> sysEnv
    String accountKey
    String accountName
    String sasToken
    Duration tokenDuration

    AzStorageOpts(Map config, Map<String,String> env=null) {
        assert config!=null
        this.sysEnv = env==null ? new HashMap<String,String>(System.getenv()) : env
        this.accountKey = config.accountKey ?: sysEnv.get('AZURE_STORAGE_ACCOUNT_KEY')
        this.accountName = config.accountName ?: sysEnv.get('AZURE_STORAGE_ACCOUNT_NAME')
        this.sasToken = config.sasToken
        this.tokenDuration = (config.tokenDuration as Duration) ?: Duration.of('12h')
    }

    Map<String,Object> getEnv() {
        Map<String, Object> props = new HashMap<>();
        props.put(AzFileSystemProvider.AZURE_STORAGE_ACCOUNT_KEY, accountKey)
        props.put(AzFileSystemProvider.AZURE_STORAGE_ACCOUNT_NAME, accountName)
        props.put(AzFileSystemProvider.AZURE_STORAGE_SAS_TOKEN, sasToken)
        return props
    }

}
