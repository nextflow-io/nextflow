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

        if (cfg.storage().accountKey && cfg.storage().sasToken) {
            throw new IllegalArgumentException("Azure Storage Access key and SAS token detected. Only one is allowed")
        }

        result.AZURE_STORAGE_ACCOUNT = cfg.storage().accountName
        // In theory, generating an impromptu SAS token for authentication methods other than
        // `azure.storage.sasToken` should not be necessary, because those methods should already allow sufficient
        // access for normal operation. Nevertheless, #5287 heavily implies that failing to do so causes the Azure
        // Storage plugin or Fusion to fail. In any case, it may be possible to remove this in the future.
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

        // For Active Directory and Managed Identity, we cannot generate an *account* SAS token, but we can generate
        // a *container* SAS token for the work directory.
        if (cfg.activeDirectory().isConfigured() || cfg.managedIdentity().isConfigured()) {
            return AzHelper.generateContainerSasWithActiveDirectory(Global.session.workDir, cfg.storage().tokenDuration)
        }

        // Shared Key authentication can use an account SAS token
        return AzHelper.generateAccountSasWithAccountKey(Global.session.workDir, cfg.storage().tokenDuration)
    }
}
