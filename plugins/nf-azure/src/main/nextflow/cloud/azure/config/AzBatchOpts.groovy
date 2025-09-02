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


import java.util.regex.Matcher
import java.util.regex.Pattern

import groovy.transform.CompileStatic
import nextflow.Global
import nextflow.Session
import nextflow.cloud.CloudTransferOptions
import nextflow.config.schema.ConfigOption
import nextflow.config.schema.ConfigScope
import nextflow.config.schema.PlaceholderName
import nextflow.fusion.FusionHelper
import nextflow.script.dsl.Description
import nextflow.util.Duration
import nextflow.util.StringUtils

/**
 * Model Azure Batch pool config settings
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class AzBatchOpts implements ConfigScope, CloudTransferOptions {

    static final private Pattern ENDPOINT_PATTERN = ~/https:\/\/(\w+)\.(\w+)\.batch\.azure\.com/

    private Map<String,String> sysEnv

    int maxParallelTransfers
    int maxTransferAttempts
    Duration delayBetweenAttempts

    @ConfigOption
    @Description("""
        The batch service account name. Defaults to environment variable `AZURE_BATCH_ACCOUNT_NAME`.
    """)
    final String accountName

    @ConfigOption
    @Description("""
        The batch service account key. Defaults to environment variable `AZURE_BATCH_ACCOUNT_KEY`.
    """)
    final String accountKey

    @ConfigOption
    @Description("""
        Enable the automatic creation of batch pools specified in the Nextflow configuration file (default: `false`).
    """)
    final Boolean allowPoolCreation

    @ConfigOption
    @Description("""
        Enable the automatic creation of batch pools depending on the pipeline resources demand (default: `true`).
    """)
    final Boolean autoPoolMode

    @ConfigOption(types=[String])
    @Description("""
        The mode in which the `azcopy` tool is installed by Nextflow (default: `'node'`).
    """)
    final CopyToolInstallMode copyToolInstallMode

    @ConfigOption
    @Description("""
        Delete all jobs when the workflow completes (default: `false`).
    """)
    final Boolean deleteJobsOnCompletion

    @ConfigOption
    @Description("""
        Delete all compute node pools when the workflow completes (default: `false`).
    """)
    final Boolean deletePoolsOnCompletion

    @ConfigOption
    @Description("""
        Delete each task when it completes (default: `true`).
    """)
    final Boolean deleteTasksOnCompletion

    @ConfigOption
    @Description("""
        The batch service endpoint e.g. `https://nfbatch1.westeurope.batch.azure.com`.
    """)
    final String endpoint

    @ConfigOption
    @Description("""
        The maximum elapsed time that jobs may run, measured from the time they are created (default: `30d`).
    """)
    final Duration jobMaxWallClockTime

    @ConfigOption
    @Description("""
        The name of the batch service region, e.g. `westeurope` or `eastus2`. Not needed when the endpoint is specified.
    """)
    final String location

    @ConfigOption
    @Description("""
        The client ID for an Azure [managed identity](https://learn.microsoft.com/en-us/entra/identity/managed-identities-azure-resources/overview) that is available on all Azure Batch node pools. This identity is used by Fusion to authenticate to Azure storage. If set to `'auto'`, Fusion will use the first available managed identity.
    """)
    final String poolIdentityClientId

    @PlaceholderName("<name>")
    final Map<String, AzPoolOpts> pools

    @ConfigOption
    @Description("""
        When the workflow completes, set all jobs to terminate on task completion (default: `true`).
    """)
    final Boolean terminateJobsOnCompletion

    @ConfigOption(types=[String])
    @Description("""
        The behaviour when Azure Batch active job quota is reached. When `'retry'` wait for job quota to clear before submitting new jobs; when `'error'` raise an error and fail (default: `'error'`).
    """)
    final JobLimitBehaviour behaviourOnJobLimit

    AzBatchOpts(Map config, Map<String,String> env=null) {
        assert config!=null
        sysEnv = env==null ? new HashMap<String,String>(System.getenv()) : env
        accountName = config.accountName ?: sysEnv.get('AZURE_BATCH_ACCOUNT_NAME')
        accountKey = config.accountKey ?: sysEnv.get('AZURE_BATCH_ACCOUNT_KEY')
        endpoint = config.endpoint
        location = config.location
        autoPoolMode = config.autoPoolMode as Boolean
        allowPoolCreation = config.allowPoolCreation as Boolean
        terminateJobsOnCompletion = config.terminateJobsOnCompletion != Boolean.FALSE
        deleteJobsOnCompletion = config.deleteJobsOnCompletion as Boolean
        deletePoolsOnCompletion = config.deletePoolsOnCompletion as Boolean
        deleteTasksOnCompletion = config.deleteTasksOnCompletion as Boolean
        jobMaxWallClockTime = config.jobMaxWallClockTime ? config.jobMaxWallClockTime as Duration : Duration.of('30d')
        poolIdentityClientId = config.poolIdentityClientId
        pools = parsePools(config.pools instanceof Map ? config.pools as Map<String,Map> : Collections.<String,Map>emptyMap())
        maxParallelTransfers = config.maxParallelTransfers ? config.maxParallelTransfers as int : MAX_TRANSFER
        maxTransferAttempts = config.maxTransferAttempts ? config.maxTransferAttempts as int : MAX_TRANSFER_ATTEMPTS
        delayBetweenAttempts = config.delayBetweenAttempts ? config.delayBetweenAttempts as Duration : DEFAULT_DELAY_BETWEEN_ATTEMPTS
        copyToolInstallMode = config.copyToolInstallMode as CopyToolInstallMode
        behaviourOnJobLimit = config.behaviourOnJobLimit as JobLimitBehaviour ?: JobLimitBehaviour.error
    }

    static Map<String,AzPoolOpts> parsePools(Map<String,Map> pools) {
        final result = new LinkedHashMap<String,AzPoolOpts>()
        for( Map.Entry<String,Map> entry : pools ) {
            result[entry.key] = new AzPoolOpts( entry.value )
        }
        if( !result.keySet().contains('auto') )
            result.put('auto', new AzPoolOpts())
        return result
    }

    AzPoolOpts pool(String name) {
        return pools.get(name)
    }

    AzPoolOpts autoPoolOpts() {
        pool('auto')
    }

    String toString() {
        "endpoint=$endpoint; account-name=$accountName; account-key=${StringUtils.redact(accountKey)}"
    }

    private List<String> endpointParts() {
        // try to infer the account name from the endpoint
        Matcher m
        if( endpoint && (m = ENDPOINT_PATTERN.matcher(endpoint)).matches() ) {
            return [ m.group(1), m.group(2) ]
        }
        else {
            return Collections.emptyList()
        }
    }

    String getAccountName() {
        if( accountName )
            return accountName
        return endpointParts()[0]
    }

    String getLocation() {
        if( location )
            return location
        // try to infer the location name from the endpoint
        return endpointParts()[1]
    }

    String getEndpoint() {
        if( endpoint )
            return endpoint
        if( accountName && location )
            return "https://${accountName}.${location}.batch.azure.com"
        return null
    }

    boolean canCreatePool() {
        allowPoolCreation || autoPoolMode
    }

    CopyToolInstallMode getCopyToolInstallMode() {
        // if the `installAzCopy` is not specified
        // `true` is returned when the pool is not create by Nextflow
        // since it can be a pool provided by the user which does not
        // provide the required `azcopy` tool
        if( copyToolInstallMode )
            return copyToolInstallMode
        if( FusionHelper.isFusionEnabled((Session) Global.session) )
            return CopyToolInstallMode.off
        canCreatePool() ? CopyToolInstallMode.node : CopyToolInstallMode.task
    }
}
