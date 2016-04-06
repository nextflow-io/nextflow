package nextflow.container

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ShifterBuilderTest extends Specification {

    def 'should return docker prefix' () {

        expect:
        ShifterBuilder.normalizeImageName(null, null) == null

        ShifterBuilder.normalizeImageName('busybox', [:]) == 'docker:busybox:latest'
        ShifterBuilder.normalizeImageName('busybox:latest', [:]) == 'docker:busybox:latest'
        ShifterBuilder.normalizeImageName('busybox:1.0', [:]) == 'docker:busybox:1.0'
        ShifterBuilder.normalizeImageName('docker:busybox:1.0', [:]) == 'docker:busybox:1.0'

        ShifterBuilder.normalizeImageName('hello/world', [:]) == 'docker:hello/world:latest'
        ShifterBuilder.normalizeImageName('hello/world:tag-name', [:]) == 'docker:hello/world:tag-name'
        ShifterBuilder.normalizeImageName('docker:hello/world:tag-name', [:]) == 'docker:hello/world:tag-name'

        ShifterBuilder.normalizeImageName('docker:busybox', [:]) == 'docker:busybox:latest'
    }


    def 'should build the shifter run command' () {

        expect:
        new ShifterBuilder('busybox').build() == 'shifter --image busybox'
        new ShifterBuilder('busybox').params(verbose: true).build() == 'shifter --verbose --image busybox'
        new ShifterBuilder('ubuntu:latest').params(entry: '/bin/bash').build() == 'shifter --image ubuntu:latest /bin/bash'

    }

}
