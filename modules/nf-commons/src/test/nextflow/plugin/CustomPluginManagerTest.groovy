package nextflow.plugin

import nextflow.Const
import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CustomPluginManagerTest extends Specification {

    def 'should NF version' () {
        given:
        def manager = Spy(CustomPluginManager)
        expect:
        manager.getSystemVersion() == Const.APP_VER
    }

    def 'should create ver manager' () {
        given:
        def manager = Spy(CustomPluginManager)
        expect:
        manager.createVersionManager() instanceof CustomVersionManager
    }
}
