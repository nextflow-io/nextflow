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

import groovy.util.logging.Slf4j
import nextflow.Global
import nextflow.cloud.azure.batch.AzHelper
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
@Slf4j
class AzFusionEnv implements FusionEnv {

    @Override
    Map<String, String> getEnvironment(String scheme, FusionConfig config) {
        if (scheme != 'az') {
            return Collections.<String, String> emptyMap()
        }

        final cfg = AzConfig.config
        final result = new LinkedHashMap(10)

        if (!cfg.storage().accountName) {
            throw new IllegalArgumentException("Missing Azure Storage account name")
        }

        result.AZURE_STORAGE_ACCOUNT = cfg.storage().accountName

        if (cfg.activeDirectory().isConfigured()) {
            result.AZURE_TENANT_ID = cfg.activeDirectory().tenantId
            result.AZURE_CLIENT_ID = cfg.activeDirectory().servicePrincipalId
            result.AZURE_CLIENT_SECRET = cfg.activeDirectory().servicePrincipalSecret
            return result
        }

        if (cfg.managedIdentity().isConfigured()) {
            // User-assigned Managed Identity
            if (cfg.managedIdentity().clientId) {
                result.AZURE_CLIENT_ID = cfg.managedIdentity().clientId
            }
            // System Managed Identity
            return result
        }

        // Shared Key authentication or Account SAS token
        result.AZURE_STORAGE_SAS_TOKEN = getOrCreateSasToken()
        return result
    }

    /**
     * Return the SAS token if it is defined in the configuration, otherwise generate one based on the requested
     * authentication method.
     */
    synchronized String getOrCreateSasToken() {

        final cfg = AzConfig.config

        // If a SAS token is already defined in the configuration, just return it
        if (cfg.storage().sasToken) {
            return cfg.storage().sasToken
        }

        // Shared Key authentication can use an account SAS token
        return AzHelper.generateAccountSasWithAccountKey(Global.session.workDir, cfg.storage().tokenDuration)
    }
}
