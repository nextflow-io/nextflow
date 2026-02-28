/*
 * Copyright 2013-2026, Seqera Labs
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

    /**
     * Get the environment variables for Fusion with Azure, taking into
     * account a managed identity ID if provided
     *
     * @param scheme The storage scheme ('az' for Azure)
     * @param config The Fusion configuration
     * @param managedIdentityId Optional managed identity client ID from pool options
     * @return Map of environment variables
     */
    @Override
    Map<String, String> getEnvironment(String scheme, FusionConfig config) {
        if (scheme != 'az') {
            return Collections.<String, String> emptyMap()
        }

        final cfg = AzConfig.config
        final managedIdentityId = cfg.batch().poolIdentityClientId
        final result = new LinkedHashMap(10)

        if (!cfg.storage().accountName) {
            throw new IllegalArgumentException("Missing Azure Storage account name")
        }

        result.AZURE_STORAGE_ACCOUNT = cfg.storage().accountName

        // If pool has a managed identity, ONLY add the MSI client ID
        // DO NOT add any SAS token or reference cfg.storage().sasToken
        if (managedIdentityId) {
            // Fusion will try and pick up a managed identity that is available.
            // We recommend explicitly setting the config item to the managed ID so you know which one is being used.
            // However if set to 'true' it will use whichever is available.
            // This can be helpful if the pools have different managed identities.
            if (managedIdentityId != 'auto') {
                result.FUSION_AZ_MSI_CLIENT_ID = managedIdentityId
            }
            // No SAS token is added or generated
            return result
        }

        // If no managed identity, use the standard environment with SAS token
        result.AZURE_STORAGE_SAS_TOKEN = getOrCreateSasToken()

        return result
    }

    /**
     * Return the SAS token if it is defined in the configuration, otherwise generate one based on the requested
     * authentication method.
     */
    synchronized String getOrCreateSasToken() {
        final cfg = AzConfig.config

        // Check for incompatible configuration
        if (cfg.storage().accountKey && cfg.storage().sasToken) {
            throw new IllegalArgumentException("Azure Storage Access key and SAS token detected. Only one is allowed")
        }

        // If a SAS token is already defined in the configuration, just return it
        if (cfg.storage().sasToken) {
            return cfg.storage().sasToken
        }

        // For Active Directory and Managed Identity, we cannot generate an *account* SAS token, but we can generate
        // a *container* SAS token for the work directory.
        if (cfg.activeDirectory().isConfigured() || cfg.managedIdentity().isConfigured()) {
            return AzHelper.generateContainerSasWithActiveDirectory(Global.session.workDir, cfg.storage().tokenDuration)
        }

        // Shared Key authentication can use an account SAS token
        return AzHelper.generateAccountSasWithAccountKey(Global.session.workDir, cfg.storage().tokenDuration)
    }
}
