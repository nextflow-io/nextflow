package nextflow.cloud.azure.batch

import nextflow.Session
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.config.AzBatchOpts
import nextflow.cloud.azure.config.AzPoolOpts
import nextflow.exception.AbortOperationException
import spock.lang.Specification

/**
 * Test for AzBatchExecutor validation logic
 */
class AzBatchExecutorTest extends Specification {

    def 'should validate low priority VMs for BatchService allocation mode'() {
        given:
        def CONFIG = [
            batch: [
                endpoint: 'https://testaccount.eastus.batch.azure.com',
                pools: [
                    'pool1': [vmType: 'Standard_D2_v2', lowPriority: true],
                    'pool2': [vmType: 'Standard_D2_v2', lowPriority: false]
                ]
            ]
        ]
        
        and:
        def config = new AzConfig(CONFIG)
        def batchService = Mock(AzBatchService) {
            getPoolAllocationMode() >> 'BATCH_SERVICE'
        }
        
        and:
        def executor = new AzBatchExecutor()
        executor.config = config
        executor.batchService = batchService

        when:
        executor.validateLowPriorityVMs()

        then:
        def e = thrown(AbortOperationException)
        e.message.contains('Azure Low Priority VMs are deprecated and no longer supported for Batch Managed pool allocation mode')
        e.message.contains('pool1')
    }

    def 'should allow low priority VMs for UserSubscription allocation mode'() {
        given:
        def CONFIG = [
            batch: [
                endpoint: 'https://testaccount.eastus.batch.azure.com',
                pools: [
                    'pool1': [vmType: 'Standard_D2_v2', lowPriority: true]
                ]
            ]
        ]
        
        and:
        def config = new AzConfig(CONFIG)
        def batchService = Mock(AzBatchService) {
            getPoolAllocationMode() >> 'USER_SUBSCRIPTION'
        }
        
        and:
        def executor = new AzBatchExecutor()
        executor.config = config
        executor.batchService = batchService

        when:
        executor.validateLowPriorityVMs()

        then:
        noExceptionThrown()
    }

    def 'should handle unknown allocation mode gracefully'() {
        given:
        def CONFIG = [
            batch: [
                endpoint: 'https://testaccount.eastus.batch.azure.com',
                pools: [
                    'pool1': [vmType: 'Standard_D2_v2', lowPriority: true]
                ]
            ]
        ]
        
        and:
        def config = new AzConfig(CONFIG)
        def batchService = Mock(AzBatchService) {
            getPoolAllocationMode() >> null
        }
        
        and:
        def executor = new AzBatchExecutor()
        executor.config = config
        executor.batchService = batchService

        when:
        executor.validateLowPriorityVMs()

        then:
        noExceptionThrown()
    }

    def 'should not validate when no low priority VMs configured'() {
        given:
        def CONFIG = [
            batch: [
                endpoint: 'https://testaccount.eastus.batch.azure.com',
                pools: [
                    'pool1': [vmType: 'Standard_D2_v2', lowPriority: false],
                    'pool2': [vmType: 'Standard_D2_v2']
                ]
            ]
        ]
        
        and:
        def config = new AzConfig(CONFIG)
        def batchService = Mock(AzBatchService) {
            getPoolAllocationMode() >> 'BATCH_SERVICE'
        }
        
        and:
        def executor = new AzBatchExecutor()
        executor.config = config
        executor.batchService = batchService

        when:
        executor.validateLowPriorityVMs()

        then:
        noExceptionThrown()
    }
}