package nextflow.script

import nextflow.Session
import nextflow.executor.Executor
import nextflow.executor.ExecutorFactory
import nextflow.processor.TaskProcessor
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ProcessDefTest extends Specification {

    def 'should clone a process with a new name '() {

        given:
        def OWNER = Mock(BaseScript)
        def BODY = new BodyDef({->}, '')
        def CONFIG = new ProcessConfig(OWNER, 'foo')
        def proc = new ProcessDef(OWNER, 'foo', CONFIG, BODY)

        when:
        def copy = proc.cloneWithName('foo_alias')
        then:
        copy.getName() == 'foo_alias'
        copy.getSimpleName() == 'foo_alias'
        copy.getBaseName() == 'foo'
        copy.getOwner() == OWNER
        copy.taskBody.class == BODY.class
        !copy.taskBody.is(BODY)

        when:
        copy = proc.cloneWithName('flow1:flow2:foo')
        then:
        copy.getName() == 'flow1:flow2:foo'
        copy.getSimpleName() == 'foo'
        copy.getBaseName() == 'foo'
        copy.getOwner() == OWNER
        copy.taskBody.class == BODY.class
        !copy.taskBody.is(BODY)
    }

    def 'should apply process config' () {
        given:
        def OWNER = Mock(BaseScript)
        def CONFIG = new ProcessConfig(OWNER, 'foo')
        def BODY = new BodyDef({->}, 'echo hello')
        def proc = new ProcessDef(OWNER, 'foo', CONFIG, BODY)
        and:
        proc.session = Mock(Session) {
            getConfig() >> [
                process: [
                    cpus: 2,
                    memory: '1GB',
                    'withName:foo': [memory: '3GB'],
                    'withName:bar': [cpus: 4, memory: '4GB'],
                    'withName:flow1:flow2:flow3:bar': [memory: '8GB']
                ]
            ]
        }

        when:
        def copy = proc.clone()
        copy.initialize()
        then:
        def cfg1 = copy.processConfig.createTaskConfig()
        cfg1.getCpus()==2           // taken from the generic config
        cfg1.getMemory().giga == 3  // taken from the `foo` config

        when:
        copy = proc.cloneWithName('flow1:bar')
        copy.initialize()
        then:
        def cfg2 = copy.processConfig.createTaskConfig()
        cfg2.getCpus()==4           // taken from the `bar` config
        cfg2.getMemory().giga == 4  // taken from the `bar` config


        when:
        copy = proc.cloneWithName('flow1:flow2:flow3:bar')
        copy.initialize()
        then:
        def cfg3 = copy.processConfig.createTaskConfig()
        cfg3.getCpus()==4           // <-- taken from `withName: foo`
        cfg3.getMemory().giga == 8  // <-- taken from `withName: 'flow1:flow2:flow3:bar'`
    }

    def 'should apply config when calling initialize multiple times (idempotent)' () {
        given:
        def OWNER = Mock(BaseScript)
        def CONFIG = new ProcessConfig(OWNER, 'foo')
        def BODY = new BodyDef({->}, 'echo hello')
        def proc = new ProcessDef(OWNER, 'foo', CONFIG, BODY)
        and:
        proc.session = Mock(Session) {
            getConfig() >> [
                process: [
                    cpus: 2,
                    'withName:foo': [cpus: 8]
                ]
            ]
        }

        when: 'initialize is called multiple times'
        proc.initialize()
        proc.initialize()
        proc.initialize()
        then: 'config is applied correctly and only once'
        def cfg = proc.processConfig.createTaskConfig()
        cfg.getCpus() == 8
    }

    def 'should apply config via createTaskProcessor without explicit initialize' () {
        given:
        def OWNER = Mock(BaseScript)
        def CONFIG = new ProcessConfig(OWNER, 'foo')
        CONFIG.container = 'source-container:1.0'  // container from source
        def BODY = new BodyDef({->}, 'echo hello')
        def proc = new ProcessDef(OWNER, 'foo', CONFIG, BODY)
        and:
        def mockExecutor = Mock(Executor) { getName() >> 'local' }
        def mockExecutorFactory = Mock(ExecutorFactory) {
            getExecutor(_, _, _, _) >> mockExecutor
        }
        def mockProcessFactory = Mock(ProcessFactory) {
            newTaskProcessor(_, _, _, _) >> Mock(TaskProcessor)
        }
        proc.session = Mock(Session) {
            getConfig() >> [
                process: [
                    'withName:foo': [container: 'config-container:2.0']  // override from config
                ]
            ]
            getExecutorFactory() >> mockExecutorFactory
            newProcessFactory(_) >> mockProcessFactory
        }

        when: 'createTaskProcessor is called without explicit initialize'
        proc.createTaskProcessor()
        then: 'config settings from withName selector are applied'
        def cfg = proc.processConfig.createTaskConfig()
        cfg.getContainer() == 'config-container:2.0'  // should be overridden by config
    }
}
