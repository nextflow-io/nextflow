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

package nextflow.cloud.azure.batch

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.cloud.azure.config.AzPoolOpts
import nextflow.cloud.azure.fusion.AzFusionEnv
import nextflow.fusion.FusionAwareTask
import nextflow.fusion.FusionScriptLauncher
import nextflow.plugin.Plugins
import nextflow.fusion.FusionEnv

/**
 * Adapter class for Azure Batch Fusion integration
 * Manages the task, launcher, and pool options
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
class AzureBatchFusionAdapter {

    private FusionAwareTask task
    private FusionScriptLauncher launcher
    private AzPoolOpts poolOpts

    AzureBatchFusionAdapter(FusionAwareTask task, FusionScriptLauncher launcher, AzPoolOpts poolOpts) {
        this.task = task
        this.launcher = launcher
        this.poolOpts = poolOpts
    }
    
    /**
     * Find the AzFusionEnv extension
     * 
     * @return The AzFusionEnv instance or null if not found
     */
    private AzFusionEnv findAzFusionEnv() {
        final extensions = Plugins.getExtensions(FusionEnv)
        if (!extensions) {
            log.warn("No FusionEnv extensions found")
            return null
        }
        
        final azFusionEnv = extensions.find { it instanceof AzFusionEnv }
        if (!azFusionEnv) {
            log.warn("AzFusionEnv extension not found among available extensions: ${extensions*.getClass().name}")
            return null
        }
        
        return (AzFusionEnv)azFusionEnv
    }
    
    /**
     * Get environment variables for the task, including Azure-specific ones
     * Filters out the SAS token when managed identity is used
     *
     * @return Map of environment variables
     */
    Map<String, String> getEnvironment() {
        // Get base environment from launcher
        final baseEnv = launcher.fusionEnv()
        if (baseEnv == null) {
            log.warn("Fusion launcher returned null environment")
            return Collections.<String,String>emptyMap()
        }
        
        final result = new LinkedHashMap<String,String>(baseEnv)
        
        // If using Azure and pool has managed identity, use AzFusionEnv to get environment
        if (poolOpts?.managedIdentityId) {
            // Find AzFusionEnv extension
            final azFusionEnv = findAzFusionEnv()
            if (azFusionEnv) {
                // Set the pool options so it has access to managedIdentityId
                azFusionEnv.setPoolOpts(poolOpts)
                
                // Get environment with managed identity
                final env = azFusionEnv.getEnvironment('az', null)
                if (env) {
                    result.putAll(env)
                }
                
                log.debug("Using managed identity for Azure Batch Fusion task: ${poolOpts.managedIdentityId}")
                
                // Make sure SAS token is removed if somehow present
                result.remove('AZURE_STORAGE_SAS_TOKEN')
            }
            else {
                // Fallback to direct environment variable if AzFusionEnv not found
                log.warn("AzFusionEnv extension not found, using direct MSI environment variable")
                result.put('FUSION_AZ_MSI_CLIENT_ID', poolOpts.managedIdentityId)
                result.remove('AZURE_STORAGE_SAS_TOKEN')
                
                // Add Azure storage account name if not already present
                if (!result.containsKey('AZURE_STORAGE_ACCOUNT')) {
                    final cfg = nextflow.cloud.azure.config.AzConfig.config
                    if (cfg.storage().accountName) {
                        result.put('AZURE_STORAGE_ACCOUNT', cfg.storage().accountName)
                    }
                }
            }
        }
        
        return result
    }
    
    /**
     * Get the command line for submitting the Fusion task
     *
     * @return List of command arguments
     */
    List<String> fusionSubmitCli() {
        return task.fusionSubmitCli()
    }
} 