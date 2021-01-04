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

}
