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

package nextflow.cloud.azure.fusion

import groovy.transform.CompileStatic
import nextflow.cloud.azure.config.AzConfig
import nextflow.fusion.FusionConfig
import nextflow.fusion.FusionEnv

/**
 * Implement environment provider for Azure specific variables
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzFusionEnv implements FusionEnv {

    @Override
    Map<String, String> getEnvironment(String scheme, FusionConfig config) {
        if( scheme!='az' )
            return Collections.<String,String>emptyMap()

        final cfg = AzConfig.config.storage()

        Map.of('AZURE_STORAGE_SAS_TOKEN', cfg.sasToken,
                'AZURE_STORAGE_ACCOUNT', cfg.accountName)
    }
}
