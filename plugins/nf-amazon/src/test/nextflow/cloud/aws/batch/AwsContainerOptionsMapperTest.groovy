package nextflow.cloud.aws.batch

import nextflow.processor.TaskConfig
import nextflow.processor.TaskRun
import spock.lang.Specification

/**
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
class AwsContainerOptionsMapperTest extends Specification {

    def 'should set env vars'() {

        given:
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(containerOptions: '--env VAR_FOO -e VAR_FOO2=value2 --env VAR_FOO3=value3')
        def handler = Spy(AwsBatchTaskHandler)
        handler.task >> task

        when:
        def job = handler.makeJobDefRequest('repo/any_image:latest')
        then:
        1 * handler.getAwsOptions() >> { new AwsOptions(cliPath: '/bin/aws') }

        def environment = job.getContainerProperties().getEnvironment()
        environment.size() == 3
        environment.get(0).toString() == '{Name: VAR_FOO,}'
        environment.get(1).toString() == '{Name: VAR_FOO3,Value: value3}'
        environment.get(2).toString() == '{Name: VAR_FOO2,Value: value2}'
    }

    def 'should set ulimits'() {

        given:
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(containerOptions: '--ulimit nofile=1280:2560 --ulimit nproc=16:32')
        def handler = Spy(AwsBatchTaskHandler)
        handler.task >> task

        when:
        def job = handler.makeJobDefRequest('repo/any_image:latest')
        then:
        1 * handler.getAwsOptions() >> { new AwsOptions(cliPath: '/bin/aws') }

        def properties = job.getContainerProperties()
        properties.getUlimits().size() == 2
        properties.getUlimits().get(0).toString() == '{HardLimit: 2560,Name: nofile,SoftLimit: 1280}'
        properties.getUlimits().get(1).toString() == '{HardLimit: 32,Name: nproc,SoftLimit: 16}'
    }

    def 'should set privileged, readonly and user'() {

        given:
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(containerOptions: '--privileged --read-only --user nf-user')
        def handler = Spy(AwsBatchTaskHandler)
        handler.task >> task

        when:
        def job = handler.makeJobDefRequest('repo/any_image:latest')
        then:
        1 * handler.getAwsOptions() >> { new AwsOptions(cliPath: '/bin/aws') }

        job.getContainerProperties().getPrivileged() == true
        job.getContainerProperties().getReadonlyRootFilesystem() == true
        job.getContainerProperties().getUser() == 'nf-user'
    }

    def 'should set tmpfs linux params'() {

        given:
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(containerOptions: '--tmpfs /run:rw,noexec,nosuid,size=64 --tmpfs /app:ro,size=128')
        def handler = Spy(AwsBatchTaskHandler)
        handler.task >> task

        when:
        def job = handler.makeJobDefRequest('repo/any_image:latest')
        then:
        1 * handler.getAwsOptions() >> { new AwsOptions(cliPath: '/bin/aws') }

        def params = job.getContainerProperties().getLinuxParameters()
        params.getTmpfs().size() == 2
        params.getTmpfs().get(0).toString() == '{ContainerPath: /run,Size: 64,MountOptions: [rw, noexec, nosuid]}'
        params.getTmpfs().get(1).toString() == '{ContainerPath: /app,Size: 128,MountOptions: [ro]}'
    }

    def 'should set memory linux params'() {

        given:
        def task = Mock(TaskRun)
        task.getName() >> 'batch-task'
        task.getConfig() >> new TaskConfig(containerOptions: '--memory-swappiness 90 --memory-swap 2048 --shm-size 1024 ')
        def handler = Spy(AwsBatchTaskHandler)
        handler.task >> task

        when:
        def job = handler.makeJobDefRequest('repo/any_image:latest')
        then:
        1 * handler.getAwsOptions() >> { new AwsOptions(cliPath: '/bin/aws') }

        def params = job.getContainerProperties().getLinuxParameters()
        params.getSharedMemorySize() == 1024
        params.getMaxSwap() == 2048
        params.getSwappiness() == 90

    }
}
