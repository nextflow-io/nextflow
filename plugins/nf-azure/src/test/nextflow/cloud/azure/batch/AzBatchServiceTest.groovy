package nextflow.cloud.azure.batch

import java.time.Instant
import java.time.temporal.ChronoUnit
import java.util.function.Predicate

import com.azure.compute.batch.models.BatchPool
import com.azure.compute.batch.models.ElevationLevel
import com.azure.identity.ManagedIdentityCredential
import com.google.common.hash.HashCode
import nextflow.Global
import nextflow.Session
import nextflow.SysEnv
import nextflow.cloud.azure.config.AzConfig
import nextflow.cloud.azure.config.AzManagedIdentityOpts
import nextflow.cloud.azure.config.AzPoolOpts
import nextflow.cloud.azure.config.AzStartTaskOpts
import nextflow.file.FileSystemPathFactory
import nextflow.processor.TaskBean
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

    def setup() {
        SysEnv.push([:])  // <-- clear the system host env
    }

    def cleanup() {
        SysEnv.pop()      // <-- restore the system host env
    }

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
        svc.computeScore(CPUS, MemoryUnit.of(MEM), MemoryUnit.of(DISK), VM) == EXPECTED

        where:
        CPUS    | MEM     | DISK    | VM                                                                 | EXPECTED
        1       | '10 MB' | '10 MB' | [numberOfCores: 1, memoryInMB: 10, resourceDiskSizeInMB: 10]       | 0.75     // Perfect match plus name adjustment (0.75)
        2       | '10 MB' | '10 MB' | [numberOfCores: 1, memoryInMB: 10, resourceDiskSizeInMB: 10]       | null     // Too many CPUs requested
        1       | '10 GB' | '10 MB' | [numberOfCores: 1, memoryInMB: 10, resourceDiskSizeInMB: 10]       | null     // Too much memory requested
        1       | '10 MB' | '10 MB' | [numberOfCores: 2, memoryInMB: 10, resourceDiskSizeInMB: 10]       | 10.75    // VM has 1 extra CPU (weighted *10)
        1       | '5 MB'  | '10 MB' | [numberOfCores: 1, memoryInMB: 10, resourceDiskSizeInMB: 10]       | 0.755    // VM has 5MB more memory plus name adjustment (0.75)
        1       | '10 MB' | '5 MB'  | [numberOfCores: 1, memoryInMB: 10, resourceDiskSizeInMB: 10]       | 0.75     // Disk difference negligible after /100
        4       | '7 GB'  | '120 GB'| [numberOfCores: 4, memoryInMB: 7168, resourceDiskSizeInMB: 122880] | 0.75     // Basic_A3 match plus name adjustment (0.75)
        4       | '7 GB'  | '120 GB'| [numberOfCores: 8, memoryInMB: 14336, resourceDiskSizeInMB: 382976]| 50.29    // Standard_A6 (worse match due to CPU difference)
    }

    def 'should find best match for northeurope' () {
        given:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        def svc = new AzBatchService(exec)
        
        when:
        def ret = svc.findBestVm('northeurope', 4, MemoryUnit.of(7168), MemoryUnit.of(122880), null)
        then:
        ret.name == 'Basic_A3'

        when:
        ret = svc.findBestVm('northeurope', 4, MemoryUnit.of(7168), MemoryUnit.of(291840),'standard_a?')
        then:
        ret.name == 'Standard_A3'

        when:
        ret = svc.findBestVm('northeurope', 4, null, MemoryUnit.of(291840), 'standard_a?')
        then:
        ret.name == 'Standard_A3'

        when:
        ret = svc.findBestVm('northeurope', 4, MemoryUnit.of(7168), MemoryUnit.of(291840), 'standard_a2,standard_a*')
        then:
        ret.name == 'Standard_A4_v2'
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
        svc.computeSlots(CPUS, MemoryUnit.of(MEM*_1GB), MemoryUnit.of(DISK*_1GB), VM_CPUS, MemoryUnit.of(VM_MEM*_1GB), MemoryUnit.of(VM_DISK*_1GB)) == EXPECTED
        where:
        CPUS  | MEM   | DISK | VM_CPUS   | VM_MEM    | VM_DISK |   EXPECTED
        1     | 1     | 1    | 1         | 1         | 1       |   1
        and:
        1     | 1     | 1    | 4         | 64        | 4       |   1         // Task needs 1 slot, 1 CPU
        1     | 8     | 1    | 4         | 64        | 4       |   1         // Task needs 1 because the memory fits into 1 cpu share
        1     | 16    | 1    | 4         | 64        | 4       |   1         // Task needs 1 because the memory fits into 1 cpu share
        1     | 32    | 1    | 4         | 64        | 4       |   2         // Task needs 2 slots because the memory is half the VM memory and requires half the CPUs
        1     | 1     | 4    | 4         | 64        | 4       |   4         // Task needs entire disk so needs 4 CPUs/slots
        1     | 1     | 4    | 4         | 64        | 8       |   2         // Task needs half the disk so needs 2 CPUs/slots
        2     | 1     | 1    | 4         | 64        | 4       |   2         // Task needs 2 CPUs
        4     | 1     | 1    | 4         | 64        | 4       |   4         // Task needs all CPUs
        2     | 64    | 1    | 4         | 64        | 4       |   4         // Task needs all memory
        1     | 1     | 1    | 256       | 256       | 256     |   1         // Max slots is 256
    }


    def 'should configure default startTask' () {
        given:
        def CONFIG = [batch:[copyToolInstallMode: 'node']]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        def svc = new AzBatchService(exec)

        when:
        def configuredStartTask = svc.createStartTask( new AzStartTaskOpts() )
        then:
        configuredStartTask.commandLine == 'bash -c "chmod +x azcopy && mkdir $AZ_BATCH_NODE_SHARED_DIR/bin/ && cp azcopy $AZ_BATCH_NODE_SHARED_DIR/bin/"'
        configuredStartTask.resourceFiles.size()==1
        configuredStartTask.resourceFiles.first().filePath == 'azcopy'
    }

    def 'should configure custom startTask' () {
        given:
        def CONFIG = [batch:[copyToolInstallMode: 'node']]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        def svc = new AzBatchService(exec)

        when:
        def configuredStartTask = svc.createStartTask( new AzStartTaskOpts(script: 'echo hello-world') )
        then:
        configuredStartTask.commandLine == 'bash -c "chmod +x azcopy && mkdir $AZ_BATCH_NODE_SHARED_DIR/bin/ && cp azcopy $AZ_BATCH_NODE_SHARED_DIR/bin/"; bash -c \'echo hello-world\''
        and:
        configuredStartTask.resourceFiles.size()==1
        configuredStartTask.resourceFiles.first().filePath == 'azcopy'
    }

    def 'should configure not install AzCopy because copyToolInstallMode is off' () {
        given:
        def CONFIG = [batch:[copyToolInstallMode: 'off']]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        def svc = new AzBatchService(exec)

        when:
        def configuredStartTask = svc.createStartTask( new AzStartTaskOpts(script: 'echo hello-world') )
        then:
        configuredStartTask.commandLine == "bash -c 'echo hello-world'"
        configuredStartTask.resourceFiles.isEmpty()
    }

    def 'should configure not install AzCopy because copyToolInstallMode is task and quote command' () {
        given:
        def CONFIG = [batch:[copyToolInstallMode: 'task']]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        def svc = new AzBatchService(exec)

        when:
        def configuredStartTask = svc.createStartTask( new AzStartTaskOpts(script: "echo 'hello-world'") )
        then:
        configuredStartTask.commandLine == "bash -c 'echo ''hello-world'''"
        configuredStartTask.resourceFiles.isEmpty()
    }

    def 'should create null startTask because no options are enabled' () {
        given:
        def CONFIG = [batch:[copyToolInstallMode: 'off']]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        def svc = new AzBatchService(exec)

        when:
        def configuredStartTask = svc.createStartTask( new AzStartTaskOpts() )
        then:
        configuredStartTask == null
    }

    def 'should configure privileged startTask' () {
        given:
        def CONFIG = [batch:[copyToolInstallMode: 'node']]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        def svc = new AzBatchService(exec)
        and:

        when:
        def configuredStartTask = svc.createStartTask( new AzStartTaskOpts(privileged: true) )
        then:
        configuredStartTask.userIdentity.autoUser.elevationLevel == ElevationLevel.ADMIN
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

    def 'should check scaling formula for low-priority' () {
        given:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        def svc = new AzBatchService(exec)

        when:
        def formula = svc.scaleFormula( new AzPoolOpts(lowPriority: true, vmCount: 3, maxVmCount: 10, scaleInterval: Duration.of('5 min')) )
        then:
        formula.contains 'interval = TimeInterval_Minute * 5;'
        formula.contains '$TargetLowPriorityNodes = lifespan < interval ? 3 : targetPoolSize;'
    }

    def 'should  check formula vars' () {
        given:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        def svc = new AzBatchService(exec)
        and:
        def opts = new AzPoolOpts(vmCount: 3, maxVmCount: 10, scaleInterval: Duration.of('5 min'))
        def now = Instant.now()

        when:
        def vars = svc.poolCreationBindings(opts, now)
        then:
        vars == [scaleInterval: 5,
                 maxVmCount: 10,
                 vmCount: 3,
                 poolCreationTime: now.truncatedTo(ChronoUnit.MICROS).toString() ]
    }

    def 'should guess vm' () {
        given:
        def LOC = 'europe'
        def TYPE = Mock(AzVmType)
        and:
        def exec = Mock(AzBatchExecutor) { getConfig() >> new AzConfig([:]) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs: [exec])

        when:
        def ret = svc.guessBestVm(LOC, 1, null, null, 'xyz')
        then:
        1 * svc.findBestVm(LOC, 1, null, null, 'xyz')  >> TYPE
        and:
        ret == TYPE

        when:
        ret = svc.guessBestVm(LOC, 8, null, null, 'xyz_*')
        then:
        1 * svc.findBestVm(LOC, 16, null, null, 'xyz_*')  >> null
        1 * svc.findBestVm(LOC, 8, null, null, 'xyz_*')  >> TYPE
        and:
        ret == TYPE

        when:
        ret = svc.guessBestVm(LOC, 8, null, null, 'xyz_?')
        then:
        1 * svc.findBestVm(LOC, 16, null, null, 'xyz_?')  >> null
        1 * svc.findBestVm(LOC, 8, null, null, 'xyz_?')  >> TYPE
        and:
        ret == TYPE

        when:
        ret = svc.guessBestVm(LOC, 16, null, null, 'xyz*')
        then:
        1 * svc.findBestVm(LOC, 16, null, null, 'xyz*')  >> TYPE
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

    def 'should create spec for autopool' () {
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
                getResourceLabels() >> [foo: 'bar']
            }
        }

        when:
        def spec = svc.specFromAutoPool(TASK)
        then:
        1 * svc.guessBestVm(LOC, CPUS, MEM, null, TYPE) >> VM
        and:
        spec.poolId == 'nf-pool-e3331cce25aa1563d6046b3de9ec2d93-Standard_X1'
        spec.metadata == [foo: 'bar']

    }

    def 'should set jobs to automatically terminate by default' () {
        given:
        def CONFIG = [:]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])
        when:
        svc.close()
        then:
        1 * svc.terminateJobs() >> null
    }

    def 'should not cleanup jobs by default' () {
        given:
        def CONFIG = [:]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])
        when:
        svc.close()
        then:
        0 * svc.cleanupJobs() >> null
    }

    def 'should cleanup jobs if specified' () {
        given:
        def CONFIG = [batch:[deleteJobsOnCompletion: true]]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(AzBatchService, constructorArgs:[exec])
        when:
        svc.close()
        then:
        1 * svc.cleanupJobs() >> null
    }

    def 'should not cleanup pools by default' () {
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

    def 'should cleanup pools with allowPoolCreation' () {
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
        result.opts
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
        1 * svc.getPool(_) >> new BatchPool(vmSize: 'Standard_D2_v2')
        and:        
        result.vmType.name == 'Standard_D2_v2'
        result.vmType.numberOfCores == 2
        and:
        result.opts.vmType == 'Standard_D2_v2'
    }

    def 'should create retry policy' () {
        given:
        def retryCfg = [delay: '100ms', maxDelay: '200ms', maxAttempts: 300]
        def CONFIG = [batch:[location: 'northeurope'], retryPolicy: retryCfg]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(new AzBatchService(exec))

        when:
        final policy = svc.retryPolicy(Mock(Predicate))
        then:
        policy.config.delay.toMillis() == 100
        policy.config.maxDelay.toMillis() == 200
        policy.config.maxAttempts == 300

    }

    def 'should create apply policy' () {
        given:
        def CONFIG = [batch:[location: 'northeurope']]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService svc = Spy(new AzBatchService(exec))

        expect:
        svc.apply(() -> 'Hello') == 'Hello'

    }

    def 'should create task for submit' () {
        given:
        Global.session = Mock(Session) { getConfig()>>[:] }
        and:
        def POOL_ID = 'my-pool'
        def SAS = '123'
        def CONFIG = [storage: [sasToken: SAS]]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService azure = Spy(new AzBatchService(exec))
        and:
        def TASK = Mock(TaskRun) {
            getHash() >> HashCode.fromInt(1)
            getContainer() >> 'ubuntu:latest'
            getConfig() >> Mock(TaskConfig)
        }
        and:
        def SPEC = new AzVmPoolSpec(poolId: POOL_ID, vmType: Mock(AzVmType), opts: new AzPoolOpts([:]))

        when:
        def result = azure.createTask(POOL_ID, 'salmon', TASK)
        then:
        1 * azure.getPoolSpec(POOL_ID) >> SPEC
        1 * azure.computeSlots(TASK, SPEC) >> 4
        1 * azure.resourceFileUrls(TASK, SAS) >> []
        1 * azure.outputFileUrls(TASK, SAS) >> []
        and:
        result.id == 'nf-01000000'
        result.requiredSlots == 4
        and:
        result.commandLine == "sh -c 'bash .command.run 2>&1 | tee .command.log'"
        and:
        result.containerSettings.imageName == 'ubuntu:latest'
        result.containerSettings.containerRunOptions == '-v /etc/ssl/certs:/etc/ssl/certs:ro -v /etc/pki:/etc/pki:ro '
    }

    def 'should create task for submit with cpu and memory' () {
        given:
        def POOL_ID = 'my-pool'
        def SAS = '123'

        def CONFIG = [storage: [sasToken: SAS]]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService azure = Spy(new AzBatchService(exec))
        def session = Mock(Session) {
            getConfig() >>[fusion:[enabled:false]]
            statsEnabled >> true
        }
        Global.session = session
        and:
        def TASK = Mock(TaskRun) {
            getHash() >> HashCode.fromInt(2)
            getContainer() >> 'ubuntu:latest'
            getConfig() >> Mock(TaskConfig) {
                getTime() >> Duration.of('24 h')
                getCpus() >> 4
                getMemory() >> MemoryUnit.of('8 GB')
            }

        }
        and:
        def SPEC = new AzVmPoolSpec(poolId: POOL_ID, vmType: Mock(AzVmType), opts: new AzPoolOpts([:]))

        when:
        def result = azure.createTask(POOL_ID, 'salmon', TASK)
        then:
        1 * azure.getPoolSpec(POOL_ID) >> SPEC
        1 * azure.computeSlots(TASK, SPEC) >> 4
        1 * azure.resourceFileUrls(TASK, SAS) >> []
        1 * azure.outputFileUrls(TASK, SAS) >> []
        and:
        result.id == 'nf-02000000'
        result.requiredSlots == 4
        and:
        result.commandLine == "sh -c 'bash .command.run 2>&1 | tee .command.log'"
        and:
        result.containerSettings.imageName == 'ubuntu:latest'
        result.containerSettings.containerRunOptions == '--cpu-shares 4096 --memory 8192m -v /etc/ssl/certs:/etc/ssl/certs:ro -v /etc/pki:/etc/pki:ro '
        and:
        Duration.of(result.constraints.maxWallClockTime.toMillis()) == TASK.config.time
    }

    def 'should create task for submit with extra options' () {
        given:
        def POOL_ID = 'my-pool'
        def SAS = '123'

        def CONFIG = [storage: [sasToken: SAS, fileShares: [file1: [mountOptions: 'mountOptions1', mountPath: 'mountPath1']]]]
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService azure = Spy(new AzBatchService(exec))
        def session = Mock(Session) {
            getConfig() >>[fusion:[enabled:false]]
            statsEnabled >> true
        }
        Global.session = session
        and:
        def TASK = Mock(TaskRun) {
            getHash() >> HashCode.fromInt(2)
            getContainer() >> 'ubuntu:latest'
            getConfig() >> Mock(TaskConfig) {
                getContainerOptions() >> '-v /foo:/foo'
                getTime() >> Duration.of('24 h')
            }

        }
        and:
        def SPEC = new AzVmPoolSpec(poolId: POOL_ID, vmType: Mock(AzVmType), opts: new AzPoolOpts([:]))

        when:
        def result = azure.createTask(POOL_ID, 'salmon', TASK)
        then:
        1 * azure.getPoolSpec(POOL_ID) >> SPEC
        1 * azure.computeSlots(TASK, SPEC) >> 4
        1 * azure.resourceFileUrls(TASK, SAS) >> []
        1 * azure.outputFileUrls(TASK, SAS) >> []
        and:
        result.id == 'nf-02000000'
        result.requiredSlots == 4
        and:
        result.commandLine == "sh -c 'bash .command.run 2>&1 | tee .command.log'"
        and:
        result.containerSettings.imageName == 'ubuntu:latest'
        result.containerSettings.containerRunOptions == '-v /etc/ssl/certs:/etc/ssl/certs:ro -v /etc/pki:/etc/pki:ro -v /mnt/batch/tasks/fsmounts/file1:mountPath1:rw -v /foo:/foo '
        and:
        Duration.of(result.constraints.maxWallClockTime.toMillis()) == TASK.config.time
    }

    def 'should create task for submit with fusion' () {
        given:
        def SAS = '1234567890' * 10
        def AZURE = [storage: [sasToken: SAS, accountName: 'my-account']]
        Global.session = Mock(Session) { getConfig()>>[fusion:[enabled:true], azure: AZURE] }
        def WORKDIR = FileSystemPathFactory.parse('az://foo/work/dir')
        and:
        def POOL_ID = 'my-pool'
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(AZURE) }
        AzBatchService azure = Spy(new AzBatchService(exec))
        and:
        def TASK = Mock(TaskRun) {
            getHash() >> HashCode.fromInt(1)
            getContainer() >> 'ubuntu:latest'
            getConfig() >> Mock(TaskConfig)
            getWorkDir() >> WORKDIR
            toTaskBean() >> Mock(TaskBean) {
                getWorkDir() >> WORKDIR
                getInputFiles() >> [:]
            }
        }
        and:
        def SPEC = new AzVmPoolSpec(poolId: POOL_ID, vmType: Mock(AzVmType), opts: new AzPoolOpts([:]))

        when:
        def result = azure.createTask(POOL_ID, 'salmon', TASK)
        then:
        1 * azure.getPoolSpec(POOL_ID) >> SPEC
        1 * azure.computeSlots(TASK, SPEC) >> 1
        1 * azure.resourceFileUrls(TASK, SAS) >> []
        1 * azure.outputFileUrls(TASK, SAS) >> []
        and:
        result.id == 'nf-01000000'
        result.requiredSlots == 1
        and:
        result.commandLine == "/usr/bin/fusion bash /fusion/az/foo/work/dir/.command.run"
        and:
        result.containerSettings.imageName == 'ubuntu:latest'
        result.containerSettings.containerRunOptions == '-v /etc/ssl/certs:/etc/ssl/certs:ro -v /etc/pki:/etc/pki:ro --privileged -e FUSION_WORK=/fusion/az/foo/work/dir -e FUSION_TAGS=[.command.*|.exitcode|.fusion.*](nextflow.io/metadata=true),[*](nextflow.io/temporary=true) -e AZURE_STORAGE_ACCOUNT=my-account -e AZURE_STORAGE_SAS_TOKEN=1234567890123456789012345678901234567890123456789012345678901234567890123456789012345678901234567890 '
    }

    @Unroll
    def 'should create user-assigned managed identity credentials token' () {
        given:
        def config = Mock(AzConfig)
        def exec = Mock(AzBatchExecutor) {getConfig() >> new AzConfig(CONFIG) }
        AzBatchService service = Spy(new AzBatchService(exec))

        when:
        def token = service.createBatchCredentialsWithManagedIdentity()
        then:
        config.managedIdentity() >> { Mock(AzManagedIdentityOpts) }
        then:
        token instanceof ManagedIdentityCredential
        (token as ManagedIdentityCredential).clientId == EXPECTED

        where:
        CONFIG                                          | EXPECTED
        [:]                                             | null
        [managedIdentity: [clientId: 'client-123']]     | 'client-123'
    }

    def 'should cache job id' () {
        given:
        def exec = Mock(AzBatchExecutor)
        def service = Spy(new AzBatchService(exec))
        and:
        def p1 = Mock(TaskProcessor)
        def p2 = Mock(TaskProcessor)
        def t1 = Mock(TaskRun) { getProcessor()>>p1 }
        def t2 = Mock(TaskRun) { getProcessor()>>p2 }
        def t3 = Mock(TaskRun) { getProcessor()>>p2 }

        when:
        def result = service.getOrCreateJob('foo',t1)
        then:
        1 * service.createJob0('foo',t1) >> 'job1'
        and:
        result == 'job1'

        // second time is cached
        when:
        result = service.getOrCreateJob('foo',t1)
        then:
        0 * service.createJob0('foo',t1) >> null
        and:
        result == 'job1'

        // changing pool id returns a new job id
        when:
        result = service.getOrCreateJob('bar',t1)
        then:
        1 * service.createJob0('bar',t1) >> 'job2'
        and:
        result == 'job2'

        // changing process returns a new job id
        when:
        result = service.getOrCreateJob('bar',t2)
        then:
        1 * service.createJob0('bar',t2) >> 'job3'
        and:
        result == 'job3'

        // change task with the same process, return cached job id
        when:
        result = service.getOrCreateJob('bar',t3)
        then:
        0 * service.createJob0('bar',t3) >> null
        and:
        result == 'job3'
    }
}
