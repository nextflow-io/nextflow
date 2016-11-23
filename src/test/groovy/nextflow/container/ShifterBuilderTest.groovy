package nextflow.container

import spock.lang.Specification

import java.nio.file.Paths

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


    def 'test shifter env'() {

        given:
        def builder = new ShifterBuilder('x')

        expect:
        builder.makeEnv('X=1').toString() == 'X=1'
        builder.makeEnv([VAR_X:1, VAR_Y: 2]).toString() == 'VAR_X=1 VAR_Y=2'
        builder.makeEnv( Paths.get('/some/file.env') ).toString() == 'BASH_ENV="/some/file.env"'
        builder.makeEnv( new File('/some/file.env') ).toString() == 'BASH_ENV="/some/file.env"'
    }

    def 'should build the shifter run command' () {

        expect:
        new ShifterBuilder('busybox')
                .build()
                .runCommand == 'shifter --image busybox'

        new ShifterBuilder('busybox')
                .params(verbose: true)
                .build()
                .runCommand == 'shifter --verbose --image busybox'

        new ShifterBuilder('ubuntu:latest')
                .params(entry: '/bin/bash')
                .build()
                .runCommand == 'shifter --image ubuntu:latest /bin/bash'

        new ShifterBuilder('ubuntu')
                .params(entry: '/bin/bash')
                .addEnv(Paths.get("/data/env_file"))
                .build()
                .runCommand == 'BASH_ENV="/data/env_file" shifter --image ubuntu /bin/bash'

        new ShifterBuilder('fedora')
                .params(entry: '/bin/bash')
                .addEnv([VAR_X:1, VAR_Y:2])
                .addEnv("VAR_Z=3")
                .build()
                .runCommand == 'VAR_X=1 VAR_Y=2 VAR_Z=3 shifter --image fedora /bin/bash'

    }


}
