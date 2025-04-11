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

package nextflow.cloud.azure.batch

import nextflow.Global
import nextflow.Session
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.config.AzPoolOpts
import nextflow.cloud.azure.config.AzStorageOpts
import nextflow.cloud.azure.fusion.AzFusionEnv
import nextflow.cloud.types.CloudMachineInfo
import nextflow.fusion.FusionAwareTask
import nextflow.fusion.FusionScriptLauncher
import nextflow.plugin.Plugins
import nextflow.processor.TaskRun
import spock.lang.Specification
import spock.util.environment.RestoreSystemProperties

/**
 * Tests for the AzureBatchFusionAdapter class
 *
 * @author Adam Talbot <adam.talbot@seqera.io>
 */
class AzureBatchFusionAdapterTest extends Specification {

    def setup() {
        // Set up a mock Global.session and AzConfig for testing
        def session = Mock(Session) {
            getConfig() >> [
                azure: [
                    storage: [
                        accountName: 'mystorageaccount'
                    ]
                ]
            ]
        }
        Global.session = session
    }

    def cleanup() {
        // Clean up the mock session
        Global.session = null
    }

    def 'should provide managed identity environment variables when available'() {
        given:
        def STORAGE_ACCOUNT = 'mystorageaccount'
        def MSI_CLIENT_ID = 'pool-managed-identity-id'
        
        def poolOpts = Stub(AzPoolOpts) {
            getManagedIdentityId() >> MSI_CLIENT_ID
        }
        
        def task = Stub(FusionAwareTask)
        
        def launcher = Stub(FusionScriptLauncher) {
            fusionEnv() >> [
                FUSION_WORK_DIR: '/fusion/work',
                AZURE_STORAGE_ACCOUNT: STORAGE_ACCOUNT,
                AZURE_STORAGE_SAS_TOKEN: 'some-sas-token'
            ]
        }
        
        // Create a stub AzFusionEnv that doesn't verify interactions
        def azFusionEnv = Stub(AzFusionEnv) {
            getEnvironment('az', _, _) >> [
                AZURE_STORAGE_ACCOUNT: STORAGE_ACCOUNT,
                FUSION_AZ_MSI_CLIENT_ID: MSI_CLIENT_ID
            ]
        }
        
        // Create adapter with stubs - no spy here
        def adapter = new AzureBatchFusionAdapter(task, launcher, poolOpts) {
            @Override
            protected AzFusionEnv findAzFusionEnv() {
                return azFusionEnv
            }
        }
        
        when: 'We get the environment'
        def env = adapter.getEnvironment()
        
        then: 'The environment has the expected values'
        env.AZURE_STORAGE_ACCOUNT == STORAGE_ACCOUNT
        env.FUSION_AZ_MSI_CLIENT_ID == MSI_CLIENT_ID
        !env.AZURE_STORAGE_SAS_TOKEN // SAS token should be removed
        env.FUSION_WORK_DIR == '/fusion/work' // Original env values should be preserved
    }
    
    def 'should fall back to direct env vars when AzFusionEnv not found'() {
        given:
        def STORAGE_ACCOUNT = 'mystorageaccount'
        def MSI_CLIENT_ID = 'pool-managed-identity-id'
        
        def poolOpts = Mock(AzPoolOpts) {
            getManagedIdentityId() >> MSI_CLIENT_ID
        }
        
        def task = Mock(FusionAwareTask) {
            getTask() >> Mock(TaskRun) {
                getId() >> 'test-task-id'
            }
        }
        
        def launcher = Mock(FusionScriptLauncher) {
            fusionEnv() >> [
                FUSION_WORK_DIR: '/fusion/work',
                AZURE_STORAGE_ACCOUNT: STORAGE_ACCOUNT,
                AZURE_STORAGE_SAS_TOKEN: 'some-sas-token'
            ]
        }
        
        def adapter = Spy(AzureBatchFusionAdapter, constructorArgs: [task, launcher, poolOpts]) {
            findAzFusionEnv() >> null
        }

        when:
        def env = adapter.getEnvironment()

        then:
        env.AZURE_STORAGE_ACCOUNT == STORAGE_ACCOUNT
        env.FUSION_AZ_MSI_CLIENT_ID == MSI_CLIENT_ID
        !env.AZURE_STORAGE_SAS_TOKEN // SAS token should be removed
        env.FUSION_WORK_DIR == '/fusion/work' // Original env values should be preserved
    }
    
    def 'should use AzConfig when storage account not in launcher env'() {
        given:
        def STORAGE_ACCOUNT = 'mystorageaccount'
        def MSI_CLIENT_ID = 'pool-managed-identity-id'
        
        def poolOpts = Mock(AzPoolOpts) {
            getManagedIdentityId() >> MSI_CLIENT_ID
        }
        
        def task = Mock(FusionAwareTask) {
            getTask() >> Mock(TaskRun) {
                getId() >> 'test-task-id'
            }
        }
        
        def launcher = Mock(FusionScriptLauncher) {
            fusionEnv() >> [
                FUSION_WORK_DIR: '/fusion/work',
                // No AZURE_STORAGE_ACCOUNT
                AZURE_STORAGE_SAS_TOKEN: 'some-sas-token'
            ]
        }
        
        def adapter = Spy(AzureBatchFusionAdapter, constructorArgs: [task, launcher, poolOpts]) {
            findAzFusionEnv() >> null
        }

        when:
        def env = adapter.getEnvironment()

        then:
        env.AZURE_STORAGE_ACCOUNT == STORAGE_ACCOUNT
        env.FUSION_AZ_MSI_CLIENT_ID == MSI_CLIENT_ID
        !env.AZURE_STORAGE_SAS_TOKEN // SAS token should be removed
        env.FUSION_WORK_DIR == '/fusion/work' // Original env values should be preserved
    }
    
    def 'should not modify environment when no managed identity in pool opts'() {
        given:
        def STORAGE_ACCOUNT = 'mystorageaccount'
        def SAS_TOKEN = 'some-sas-token'
        
        def poolOpts = Mock(AzPoolOpts) {
            getManagedIdentityId() >> null // No managed identity
        }
        
        def task = Mock(FusionAwareTask)
        
        def launcher = Mock(FusionScriptLauncher) {
            fusionEnv() >> [
                FUSION_WORK_DIR: '/fusion/work',
                AZURE_STORAGE_ACCOUNT: STORAGE_ACCOUNT,
                AZURE_STORAGE_SAS_TOKEN: SAS_TOKEN
            ]
        }
        
        def adapter = new AzureBatchFusionAdapter(task, launcher, poolOpts)

        when:
        def env = adapter.getEnvironment()

        then:
        env.AZURE_STORAGE_ACCOUNT == STORAGE_ACCOUNT
        env.AZURE_STORAGE_SAS_TOKEN == SAS_TOKEN // SAS token should be preserved
        env.FUSION_WORK_DIR == '/fusion/work' // Original env values should be preserved
        !env.FUSION_AZ_MSI_CLIENT_ID // No MSI client ID
    }

    def 'should use SAS token when no managed identity is available'() {
        given:
        def STORAGE_ACCOUNT = 'mystorageaccount'
        def SAS_TOKEN = 'some-sas-token'
        
        def poolOpts = Stub(AzPoolOpts) {
            getManagedIdentityId() >> null // No managed identity
        }
        
        def task = Stub(FusionAwareTask)
        
        def launcher = Stub(FusionScriptLauncher) {
            fusionEnv() >> [
                FUSION_WORK_DIR: '/fusion/work',
                AZURE_STORAGE_ACCOUNT: STORAGE_ACCOUNT,
                AZURE_STORAGE_SAS_TOKEN: SAS_TOKEN
            ]
        }
        
        // Create a stub AzFusionEnv that doesn't verify interactions
        def azFusionEnv = Stub(AzFusionEnv) {
            getEnvironment('az', _, _) >> [
                AZURE_STORAGE_ACCOUNT: STORAGE_ACCOUNT,
                AZURE_STORAGE_SAS_TOKEN: SAS_TOKEN
            ]
        }
        
        // Create adapter with stubs
        def adapter = new AzureBatchFusionAdapter(task, launcher, poolOpts) {
            @Override
            protected AzFusionEnv findAzFusionEnv() {
                return azFusionEnv
            }
        }
        
        when: 'We get the environment'
        def env = adapter.getEnvironment()
        
        then: 'The environment has the expected values'
        env.AZURE_STORAGE_ACCOUNT == STORAGE_ACCOUNT
        env.AZURE_STORAGE_SAS_TOKEN == SAS_TOKEN
        env.FUSION_WORK_DIR == '/fusion/work' // Original env values should be preserved
        !env.FUSION_AZ_MSI_CLIENT_ID // No managed identity
    }
} 