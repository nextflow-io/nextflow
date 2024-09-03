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

package nextflow.cloud.azure.fusion

import groovy.transform.CompileStatic
import nextflow.cloud.azure.config.AzConfig
import nextflow.fusion.FusionConfig
import nextflow.fusion.FusionEnv
import org.pf4j.Extension
/**
 * Implement environment provider for Azure specific variables
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Extension
@CompileStatic
class AzFusionEnv implements FusionEnv {

    @Override
    Map<String, String> getEnvironment(String scheme, FusionConfig config) {
        if (scheme != 'az')
            return Collections.<String, String> emptyMap()

        final cfg = AzConfig.config
        final result = new LinkedHashMap(10)

        if (!cfg.storage().accountName)
            throw new IllegalArgumentException("Missing Azure Storage account name")

        if (cfg.storage().accountKey && cfg.storage().sasToken)
            throw new IllegalArgumentException("Azure Storage Access key and SAS token detected. Only one is allowed")

        // the account name is always required
        result.AZURE_STORAGE_ACCOUNT = cfg.storage().accountName

        // If a Managed Identity or Service Principal is configured, Fusion only needs to know the account name
        if (cfg.managedIdentity().isConfigured() || cfg.activeDirectory().isConfigured()) {
            return result
        }

        // If a SAS token is configured, instead, Fusion also requires the token value
        if (cfg.storage().sasToken) {
            result.AZURE_STORAGE_SAS_TOKEN = cfg.storage().getOrCreateSasToken()
        }

        return result
    }
}
