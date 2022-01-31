package nextflow.plugin

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PluginSpecTest extends Specification {

    def 'should parse plugins spec' () {

        when:
        def spec = PluginSpec.parse(FQID)
        then:
        spec.id == ID
        spec.version == VER

        
        where:
        FQID        | ID    | VER
        'foo@1.0'   | 'foo' | '1.0'
    }

    def 'should compare specs' () {
        given:
        def spec1 = new PluginSpec('nf-alpha','1.0.0')
        def spec2 = new PluginSpec('nf-alpha','1.0.0')
        def spec3 = new PluginSpec('nf-delta','2.0.0')

        expect:
        spec1 == spec2
        spec1 != spec3
        and:
        spec1 < spec3
        and:
        spec1.hashCode() == spec2.hashCode()
        spec1.hashCode() != spec3.hashCode()

    }

}
