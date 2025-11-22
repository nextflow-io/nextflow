package nextflow.cloud.aws.batch

import nextflow.util.CmdLineHelper
import software.amazon.awssdk.services.batch.model.Tmpfs
import software.amazon.awssdk.services.batch.model.Ulimit
import spock.lang.Specification

/**
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
class AwsContainerOptionsMapperTest extends Specification {

    def 'should set env vars'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--env VAR_FOO -e VAR_FOO2=value2 --env VAR_FOO3=value3')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        def environment = properties.environment
        environment.size() == 3
        environment.get(0).name() == 'VAR_FOO'
        environment.get(0).value() == null
        environment.get(1).name() == 'VAR_FOO3'
        environment.get(1).value() == 'value3'
        environment.get(2).name() == 'VAR_FOO2'
        environment.get(2).value() == 'value2'
    }

    def 'should set ulimits'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--ulimit nofile=1280:2560 --ulimit nproc=16:32')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.ulimits.size() == 2
        properties.ulimits.get(0) == Ulimit.builder().hardLimit(2560).name('nofile').softLimit(1280).build()
        properties.ulimits.get(1) == Ulimit.builder().hardLimit(32).name('nproc').softLimit(16).build()

    }

    def 'should set user'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--user nf-user')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.user == 'nf-user'
    }

    def 'should set privileged'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--privileged')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.privileged
    }

    def 'should set readonly'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--read-only')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.readonlyRootFilesystem
    }

    def 'should set env'() {
        when:
        def map = CmdLineHelper.parseGnuArgs('-e x=y')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.environment.get(0).name()=='x'
        properties.environment.get(0).value()=='y'
    }

    def 'should set tmpfs linux params'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--tmpfs /run:rw,noexec,nosuid,size=64 --tmpfs /app:ro,size=128')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.linuxParameters.tmpfs().get(0) == Tmpfs.builder().containerPath('/run').size(64).mountOptions(['rw', 'noexec', 'nosuid']).build()
        properties.linuxParameters.tmpfs().get(1) == Tmpfs.builder().containerPath('/app').size(128).mountOptions(['ro']).build()
    }

    def 'should set memory swap '() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--memory-swap 2048')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.linuxParameters.maxSwap() == 2048
    }

    def 'should set shared memory size'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--shm-size 12048024')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.linuxParameters.sharedMemorySize() == 11
    }

    def 'should set shared memory size with unit in MiB'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--shm-size 256m')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.linuxParameters.sharedMemorySize() == 256
    }

    def 'should set shared memory size with unit in GiB'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--shm-size 1g')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.linuxParameters.sharedMemorySize() == 1024
    }

    def 'should set memory swappiness'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--memory-swappiness 12048024')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.linuxParameters.swappiness() == 12048024
    }

    def 'should set init'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--init')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.linuxParameters.initProcessEnabled()
    }

    def 'should set no params'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.linuxParameters == null
        properties.ulimits == null
        properties.privileged == false
        properties.readonlyRootFilesystem == false
        properties.user == null
    }
}

