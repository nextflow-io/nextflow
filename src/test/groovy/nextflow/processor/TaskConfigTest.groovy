package nextflow.processor

import nextflow.util.Duration
import nextflow.util.MemoryUnit
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskConfigTest extends Specification {

    def 'test defaults' () {

        setup:
        def config = new TaskConfig()

        expect:
        config.shell ==  ['/bin/bash','-ue']
        config.cacheable
        config.validExitCodes == [0]
        config.errorStrategy == ErrorStrategy.TERMINATE
        config.inputs instanceof Map
        config.outputs instanceof Map

    }

    def 'test setting properties' () {

        setup:
        def config = new TaskConfig()

        // setting property using method without brackets
        when:
        config.hola 'val 1'
        then:
        config.hola == 'val 1'

        // setting list values
        when:
        config.hola 1,2,3
        then:
        config.hola == [1,2,3]

        // setting named parameters attribute
        when:
        config.hola field1:'val1', field2: 'val2'
        then:
        config.hola == [field1:'val1', field2: 'val2']

        // maxDuration property
        when:
        config.maxDuration '1h'
        then:
        config.maxDuration == new Duration('1h')

        // maxMemory property
        when:
        config.maxMemory '2GB'
        then:
        config.maxMemory == new MemoryUnit('2GB')

        // generic value assigned like a 'plain' property
        when:
        config.hola = 99
        then:
        config.hola == 99

    }


    def 'test NO missingPropertyException' () {

        when:
        def config = new TaskConfig()
        def x = config.hola

        then:
        x == null
        noExceptionThrown()

    }

    def 'test MissingPropertyException' () {
        when:
        def config = new TaskConfigWrapper(new TaskConfig())
        def x = config.hola

        then:
        thrown(MissingPropertyException)
    }


    def 'test check property existence' () {

        setup:
        def config = new TaskConfig()

        expect:
        config.containsKey('echo')
        config.containsKey('shell')
        config.containsKey('validExitCodes')
        config.containsKey('inputs')
        config.containsKey('outputs')
        !config.containsKey('xyz')
        !config.containsKey('maxForks')
        config.maxForks == null


    }


}
