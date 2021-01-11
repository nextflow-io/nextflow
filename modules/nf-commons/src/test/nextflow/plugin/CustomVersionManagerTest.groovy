package nextflow.plugin

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CustomVersionManagerTest extends Specification {

    @Unroll
    def 'should compare version' () {
        given:
        def manager = new CustomVersionManager()

        expect:
        manager.compareVersions(V1, V2) == RET

        where:
        V1          | V2            | RET
        '1.0.0'     | '1.0.0'       | 0
        '1.1.0'     | '1.0.0'       | 1
        '1.0.0'     | '1.0.1'       | -1
        and:
        '20.01.1'   | '20.01.1'     | 0
        '20.02.1'   | '20.01.1'     | 1
        '20.02.1'   | '20.02.2'     | -1
        and:
        '20.01.1-edge'   | '20.01.1-edge'     | 0
        '20.02.1-edge'   | '20.01.1-edge'     | 1
        '20.02.1-edge'   | '20.02.2-edge'     | -1

    }

    @Unroll
    def 'check version constraint' () {
        given:
        def manager = new CustomVersionManager()

        expect:
        manager.checkVersionConstraint(VER, REQ) == RET

        where:
        VER                 | REQ               | RET
        '20.01.0'           | '>=20.01.0'       | true
        '20.01.0'           | '>=20.02.0'       | false
        '20.01.0-SNAPSHOT'  | '>=20.01.0'       | true
        and:
        '20.01.0-edge'      | '>=20.01.0-edge'  | true
        '20.01.0-edge'      | '>=20.02.0-edge'  | false
        and:
        '20.01.0-SNAPSHOT'  | '>=20.01.0-edge'  | true
        '20.01.0-SNAPSHOT'  | '>=20.01.0'       | true
    }

}
