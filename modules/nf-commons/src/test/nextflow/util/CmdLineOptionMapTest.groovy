package nextflow.util

import spock.lang.Specification
/**
 *
 * @author Manuele Simi <manuele.simi@gmail.com>
 */
class CmdLineOptionMapTest extends Specification {

    def 'test boolean values' () {

        setup:
        CmdLineOptionMap options = new CmdLineOptionMap()
        options.addOption('iamtrue', 'true')
        options.addOption('iamfalse', 'false')

        expect:
        options.hasOptions()
        options.exists('iamtrue')
        options.exists('iamfalse')
        options.getFirstValue('iamtrue') == 'true'
        options.getFirstValue('iamfalse') == 'false'
        options.getFirstValue('idontexist') == null

    }

    def 'test multiple values' () {

        setup:
        CmdLineOptionMap options = new CmdLineOptionMap()
        options.addOption('firstkey', 'value1')
        options.addOption('firstkey', 'value2')
        options.addOption('firstkey', 'value3')
        options.addOption('secondkey', 'value4')
        options.addOption('secondkey', 'value5')

        expect:
        options.hasOptions()
        options.exists('firstkey')
        options.getFirstValue('firstkey') == 'value1'
        options.getValues('firstkey').size() == 3
        options.getValues('firstkey').get(0) == 'value1'
        options.getValues('firstkey').get(1) == 'value2'
        options.getValues('firstkey').get(2) == 'value3'
        options.exists('secondkey')
        options.getValues('secondkey').size() == 2
        options.getValues('secondkey').get(0) == 'value4'
        options.getValues('secondkey').get(1) == 'value5'

    }

    def 'test default values' () {
        setup:
        CmdLineOptionMap options = new CmdLineOptionMap()
        options.addOption('key', '')
        options.addOption('key2', null)

        expect:
        options.exists('key')
        options.exists('key2')
        options.getFirstValueOrDefault('key', 'alt' ) == 'alt'
        options.getFirstValueOrDefault('key2', 'alt2' ) == 'alt2'
        !options.emptyOption().hasOptions()

    }

    def 'should check groovy truth' () {
        expect:
        // empty => evaluates to false
        !new CmdLineOptionMap()
        and:
        // not empty => evaluates to true
        new CmdLineOptionMap().addOption('foo','bar')
    }

    def 'should validate equals and hashcode' () {
        given:
        def map1 = CmdLineOptionMap.fromMap([foo: 'hello'])
        def map2 = CmdLineOptionMap.fromMap([foo: 'hello'])
        def map3 = CmdLineOptionMap.fromMap([foo: 'world'])

        expect:
        map1 == map2 
        map1 != map3
        and:
        map1.hashCode() == map2.hashCode()
        map1.hashCode() != map3.hashCode()
    }
}
