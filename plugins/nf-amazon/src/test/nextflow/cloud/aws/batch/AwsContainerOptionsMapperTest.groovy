package nextflow.cloud.aws.batch

import nextflow.util.CmdLineHelper
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
        def environment = properties.getEnvironment()
        environment.size() == 3
        environment.get(0).toString() == '{Name: VAR_FOO,}'
        environment.get(1).toString() == '{Name: VAR_FOO3,Value: value3}'
        environment.get(2).toString() == '{Name: VAR_FOO2,Value: value2}'
    }

    def 'should set ulimits'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--ulimit nofile=1280:2560 --ulimit nproc=16:32')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getUlimits().size() == 2
        properties.getUlimits().get(0).toString() == '{HardLimit: 2560,Name: nofile,SoftLimit: 1280}'
        properties.getUlimits().get(1).toString() == '{HardLimit: 32,Name: nproc,SoftLimit: 16}'

    }

    def 'should set user'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--user nf-user')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getUser() == 'nf-user'
    }

    def 'should set privileged'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--privileged')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getPrivileged()
    }

    def 'should set readonly'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--read-only')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getReadonlyRootFilesystem()
    }

    def 'should set env'() {
        when:
        def map = CmdLineHelper.parseGnuArgs('-e x=y')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getEnvironment().get(0).getName()=='x'
        properties.getEnvironment().get(0).getValue()=='y'
    }

    def 'should set tmpfs linux params'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--tmpfs /run:rw,noexec,nosuid,size=64 --tmpfs /app:ro,size=128')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getLinuxParameters().getTmpfs().get(0).toString() == '{ContainerPath: /run,Size: 64,MountOptions: [rw, noexec, nosuid]}'
        properties.getLinuxParameters().getTmpfs().get(1).toString() == '{ContainerPath: /app,Size: 128,MountOptions: [ro]}'
    }

    def 'should set memory swap '() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--memory-swap 2048')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getLinuxParameters().getMaxSwap() == 2048
    }

    def 'should set shared memory size'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--shm-size 12048024')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getLinuxParameters().getSharedMemorySize() == 11
    }

    def 'should set shared memory size with unit in MiB'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--shm-size 256m')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getLinuxParameters().getSharedMemorySize() == 256
    }

    def 'should set shared memory size with unit in GiB'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--shm-size 1g')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getLinuxParameters().getSharedMemorySize() == 1024
    }

    def 'should set memory swappiness'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--memory-swappiness 12048024')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getLinuxParameters().getSwappiness() == 12048024
    }

    def 'should set init'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('--init')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getLinuxParameters().getInitProcessEnabled()
    }

    def 'should set no params'() {

        when:
        def map = CmdLineHelper.parseGnuArgs('')
        def properties = AwsContainerOptionsMapper.createContainerProperties(map)
        then:
        properties.getLinuxParameters() == null
        properties.getUlimits() == null
        properties.getPrivileged() == null
        properties.getReadonlyRootFilesystem() == null
        properties.getUser() == null
    }
}

