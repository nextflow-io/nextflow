package nextflow.cloud.azure.batch

import com.microsoft.azure.batch.protocol.models.CloudPool
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.config.AzPoolOpts
import nextflow.processor.TaskConfig
import nextflow.processor.TaskProcessor
import nextflow.processor.TaskRun
import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AzBatchServiceTest extends Specification {

    static long _1GB = 1024 * 1024 * 1024

    def 'should make job id'() {
        given:
        def task = Mock(TaskRun) {
            getProcessor() >> Mock(TaskProcessor) {
                getName() >> NAME
            }
        }
        and:
        def exec = Mock(AzBatchExecutor) {
            getConfig() >> new AzConfig([:])
        }
        and:
        def svc = new AzBatchService(exec)

        expect:
        svc.makeJobId(task) =~ EXPECTED

        where:
        NAME        | EXPECTED
        'foo'       | /job-\w+-foo/
        'foo  bar'  | /job-\w+-foo_bar/
    }

    def 'should list locations' () {
        given:
        def exec = Mock(AzBatchExecutor) {
            getConfig() >> new AzConfig([:])
        }
        def svc = new AzBatchService(exec)

        when:
        def list = svc.listLocationNames()
        then:
        'eastus' in list
        'northeurope' in list
    }

    def 'should list vm names for location' () {
        given:
        def exec = Mock(AzBatchExecutor) {
            getConfig() >> new AzConfig([:])
        }
        def svc = new AzBatchService(exec)

        when:
        def names = svc.listVmNames('northeurope')
        then:
        'Standard_D1_v2' in names
        'Standard_D5_v2' in names
    }

    def 'should get size for vm' () {
        given:
        def exec = Mock(AzBatchExecutor) {
            getConfig() >> new AzConfig([:])
        }
        def svc = new AzBatchService(exec)

        when:
        def vm = svc.getVmType('northeurope', 'Standard_D1_v2')
        then:
        vm.name == 'Standard_D1_v2'
        vm.numberOfCores == 1
        vm.maxDataDiskCount == 4
        vm.memory.toMega() == 3584
        vm.osDiskSize.toMega() == 1047552
        vm.resourceDiskSize.toMega() == 51200


        when:
        def vm1= svc.getVmType('northeurope', 'standard_d1_v2')
        def vm2= svc.getVmType('northeurope', 'STANDARD_D1_V2')
        then:
        vm == vm1
        vm == vm2

    }

    @Unroll
    def 'should compute vm score' () {
        given:
        def exec = Mock(AzBatchExecutor) {
            getConfig() >> new AzConfig([:])
        }
        def svc = new AzBatchService(exec)


        expect:
        svc.computeScore(CPUS, MemoryUnit.of(MEM), VM) == EXPECTED

        
        where:
        CPUS    | MEM       | VM                                            | EXPECTED
        1       | '10 MB'   | [numberOfCores: 1, memoryInMb: 10]            | 0.0
        2       | '10 MB'   | [numberOfCores: 1, memoryInMb: 10]            | null
        1       | '10 GB'   | [numberOfCores: 1, memoryInMb: 10]            | null
        1       | '10 MB'   | [numberOfCores: 2, memoryInMb: 10]            | 1.0
        1       | '10 GB'   | [numberOfCores: 1, memoryInMb: 10]            | null
        1       | '10 GB'   | [numberOfCores: 1, memoryInMb: 11 * 1024]     | 1.0
    }

    def 'should find best match for northeurope' () {
        given:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        def svc = new AzBatchService(exec)
        
        when:
        def ret = svc.findBestVm('northeurope', 4, MemoryUnit.of(7168), null)
        then:
        ret.name == 'Basic_A3'

        when:
        ret = svc.findBestVm('northeurope', 4, MemoryUnit.of(7168), 'standard_a?')
        then:
        ret.name == 'Standard_A3'

        when:
        ret = svc.findBestVm('northeurope', 4, null, 'standard_a?')
        then:
        ret.name == 'Standard_A6'

        when:
        ret = svc.findBestVm('northeurope', 4, MemoryUnit.of(7168), 'standard_a2,standard_a*')
        then:
        ret.name == 'Standard_A3'
    }

    def 'should match familty' () {
        given:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        def svc = new AzBatchService(exec)

        expect:
        svc.matchType(FAMILY, VM_TYPE) == EXPECTED

        where:
        FAMILY      | VM_TYPE       | EXPECTED
        null        | 'Basic_A3'    | true
        'Basic_A3'  | 'Basic_A3'    | true
        'BASIC_A3'  | 'Basic_A3'    | true
        'Basic_A4'  | 'Basic_A3'    | false
        'Basic'     | 'Basic_A3'    | false
        'BASIC'     | 'Basic_A3'    | false
        'Other'     | 'Basic_A3'    | false
        'basic*'    | 'Basic_A3'    | true
        'basic*'    | 'Other_A3'    | false
        'Basic_A?'  | 'Basic_A3'    | true
        'Basic_A?'  | 'Basic_A33'   | false

    }

    @Unroll
    def 'should compute mem slots' () {
        given:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        def svc = new AzBatchService(exec)

        expect:
        svc.memSlots(TASK_MEM, VM_MEM, VM_CPUS) == SLOTS

        where:
        TASK_MEM | VM_MEM    | VM_CPUS   | SLOTS
        1.0      | 1.0f      | 1         | 1
        1.1      | 1.0f      | 1         | 2
        1.9      | 1.0f      | 1         | 2
        2.0      | 1.0f      | 1         | 2
        2.0      | 2.0f      | 1         | 1
        8.0      | 4.0f      | 1         | 2
        8.0      | 8.0f      | 4         | 4
        4.0      | 8.0f      | 4         | 2
        5.0      | 8.0f      | 4         | 3
    }

    @Unroll
    def 'should compute slots' () {
        given:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        def svc = new AzBatchService(exec)
        
        expect:
        svc.computeSlots(CPUS, MemoryUnit.of(MEM * _1GB), VM_CPUS, MemoryUnit.of(VM_MEM*_1GB)) == EXPECTED

        where:
        CPUS  | MEM   | VM_CPUS   | VM_MEM    |   EXPECTED
        1     | 1     | 1         | 1         |   1
        and:
        1     | 1     | 4         | 64        |   1
        1     | 8     | 4         | 64        |   1
        1     | 16    | 4         | 64        |   1
        1     | 32    | 4         | 64        |   2
        2     | 1     | 4         | 64        |   2
        4     | 1     | 4         | 64        |   4
        2     | 64    | 4         | 64        |   4

    }


    def 'should check scaling formula' () {
        given:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        def svc = new AzBatchService(exec)

        when:
        def formula = svc.scaleFormula( new AzPoolOpts(vmCount: 3, maxVmCount: 10, scaleInterval: Duration.of('5 min')) )
        then:
        formula.contains 'interval = TimeInterval_Minute * 5;'
        formula.contains '$TargetDedicatedNodes = lifespan < interval ? 3 : targetPoolSize;'
    }

    def 'should guess vm' () {
        given:
        def LOC = 'europe'
        def TYPE = Mock(AzVmType)
        and:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs: [exec])

        when:
        def ret = svc.guessBestVm(LOC, 1, null, 'xyz')
        then:
        1 * svc.findBestVm(LOC, 1, null, 'xyz')  >> TYPE
        and:
        ret == TYPE

        when:
        ret = svc.guessBestVm(LOC, 8, null, 'xyz_*')
        then:
        1 * svc.findBestVm(LOC, 16, null, 'xyz_*')  >> null
        1 * svc.findBestVm(LOC, 8, null, 'xyz_*')  >> TYPE
        and:
        ret == TYPE

        when:
        ret = svc.guessBestVm(LOC, 8, null, 'xyz_?')
        then:
        1 * svc.findBestVm(LOC, 16, null, 'xyz_?')  >> null
        1 * svc.findBestVm(LOC, 8, null, 'xyz_?')  >> TYPE
        and:
        ret == TYPE

        when:
        ret = svc.guessBestVm(LOC, 16, null, 'xyz*')
        then:
        1 * svc.findBestVm(LOC, 16, null, 'xyz*')  >> TYPE
        and:
        ret == TYPE

    }

    def 'should check poolid' () {
        given:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs: [exec])

        when:
        svc.checkPoolId('abc')
        then:
        noExceptionThrown()

        when:
        svc.checkPoolId('abc_10')
        then:
        noExceptionThrown()

        when:
        svc.checkPoolId('abc-10')
        then:
        noExceptionThrown()

        when:
        svc.checkPoolId('abc 10')
        then:
        thrown(IllegalArgumentException)

    }

    def 'should create spec for autotask' () {
        given:
        def LOC = 'europe'
        def CFG = new AzConfig([batch: [location: LOC]])
        def CPUS = 2
        def MEM = MemoryUnit.of('1 GB')
        def TYPE = 'Standard_X1'
        def VM = new AzVmType(name: TYPE, numberOfCores: CPUS)
        and:
        def exec = Mock(AzBatchExecutor) { getConfig() >> CFG }
        AzBatchService svc = Spy(AzBatchService, constructorArgs: [exec])
        and:
        def TASK = Mock(TaskRun) {
            getConfig() >> Mock(TaskConfig) {
                getMemory() >> MEM
                getCpus() >> CPUS
                getMachineType() >> TYPE
            }
        }

        when:
        def spec = svc.specFromAutoPool(TASK)
        then:
        1 * svc.guessBestVm(LOC, CPUS, MEM, TYPE) >> VM
        and:
        spec.poolId == 'nf-pool-c1575b88d347ebc5c068777405865b39-Standard_X1'

    }


    def 'should cleanup jobs by default' () {
        given:
        def CONFIG = [:]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])
        when:
        svc.close()
        then:
        1 * svc.cleanupJobs() >> null
    }

    def 'should cleanup jobs no cleanup jobs' () {
        given:
        def CONFIG = [batch:[deleteJobsOnCompletion: false]]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])
        when:
        svc.close()
        then:
        0 * svc.cleanupJobs() >> null
    }

    def 'should cleanup not cleanup pools by default' () {
        given:
        def CONFIG = [:]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])
        when:
        svc.close()
        then:
        0 * svc.cleanupPools() >> null
    }

    def 'should cleanup pools with autoPoolMode' () {
        given:
        def CONFIG = [batch:[autoPoolMode: true, deletePoolsOnCompletion: true]]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])
        when:
        svc.close()
        then:
        1 * svc.cleanupPools() >> null
    }

    def 'should cleanup cleanup pools with allowPoolCreation' () {
        given:
        def CONFIG = [batch:[allowPoolCreation: true, deletePoolsOnCompletion: true]]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])
        when:
        svc.close()
        then:
        1 * svc.cleanupPools() >> null
    }


    def 'should not cleanup pools without autoPoolMode or allowPoolCreation' () {
        given:
        def CONFIG = [batch:[deletePoolsOnCompletion: true]]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])
        when:
        svc.close()
        then:
        0 * svc.cleanupPools() >> null
    }

    def 'should get spec from pool config' () {
        given:
        def POOL_ID = 'foo'
        def CONFIG = [batch:[location: 'northeurope', pools: [(POOL_ID): [vmType: 'Standard_D2_v2']]]]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])

        when:
        def result = svc.specFromPoolConfig(POOL_ID)
        then:
        result.vmType.name == 'Standard_D2_v2'
        result.vmType.numberOfCores == 2
        and:
        0 * svc.getPool(_) >> null
    }

    def 'should get spec from existing pool' () {
        given:
        def POOL_ID = 'foo'
        def CONFIG = [batch:[location: 'northeurope']]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])

        when:
        def result = svc.specFromPoolConfig(POOL_ID)
        then:
        1 * svc.getPool(_) >> new CloudPool(vmSize: 'Standard_D2_v2')
        and:        
        result.vmType.name == 'Standard_D2_v2'
        result.vmType.numberOfCores == 2
    }

}
