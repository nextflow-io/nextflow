package nextflow.container

import spock.lang.Specification

import java.nio.file.Paths
import java.nio.file.Files

/**
 * @author Emilio Palumbo <emiliopalumbo@gmail.com>
 */
class ContainerHandlerTest extends Specification {

    def 'test docker is absolute image name' () {

        expect:
        ContainerHandler.isAbsoluteDockerName(image) == expected

        where:
        image                      | expected
        'hello'                    | false
        'image/name'               | false
        'registry:5000/image/name' | true
        'd.reg/image/name'         | true
        'd.reg/image'              | true

    }

    def 'test normalize docker image name' () {

        given:
        def n = new ContainerHandler([registry: registry])

        expect:
        n.normalizeDockerImageName(image) == expected

        where:
        image                       | registry  | expected
        null                        | null      | null
        null                        | 'd.reg'   | null
        'hello'                     | null      | 'hello'
        'cbcrg/hello'               | null      | 'cbcrg/hello'
        'cbcrg/hello'               | 'd.reg'   | 'd.reg/cbcrg/hello'
        'cbcrg/hello'               | 'd.reg/'  | 'd.reg/cbcrg/hello'
        'registry:5000/cbcrg/hello' | 'd.reg'   | 'registry:5000/cbcrg/hello'

    }

    def 'test normalize shifter image name' () {

        given:
        def n = new ContainerHandler([:])

        expect:
        n.normalizeShifterImageName(image) == expected

        where:
        image                         | expected
        'busybox'                     | 'docker:busybox:latest'
        'busybox:latest'              | 'docker:busybox:latest'
        'busybox:1.0'                 | 'docker:busybox:1.0'
        'docker:busybox:1.0'          | 'docker:busybox:1.0'
        'hello/world'                 | 'docker:hello/world:latest'
        'hello/world:tag-name'        | 'docker:hello/world:tag-name'
        'docker:hello/world:tag-name' | 'docker:hello/world:tag-name'
        'docker:busybox'              | 'docker:busybox:latest'
    }

    def 'test normalize singularity image name' () {

        given:
        def n = new ContainerHandler([:])

        expect:
        n.normalizeSingularityImageName(image) == expected

        where:
        image                      | expected
        null                       | null
        ''                         | null
        '/abs/path/bar.img'        | '/abs/path/bar.img'
        'file:///abs/path/bar.img' | '/abs/path/bar.img'
        'docker://library/busybox' | 'docker://library/busybox'
        'shub://busybox'           | 'shub://busybox'
        'foo'                      | 'docker://foo'
        'foo:2.0'                  | 'docker://foo:2.0'
        'foo.img'                  | 'docker://foo.img'
    }

    def 'test singularity relative path exists' () {

        setup:
        def base = Files.createTempDirectory('test')
        def foo = base.resolve('foo'); foo.mkdir()
        def bar = Files.createFile(foo.resolve('bar'))
        def img = Files.createFile(base.resolve('bar.img'))
        def n = new ContainerHandler([:], base)

        expect:
        n.normalizeSingularityImageName('foo/bar') == bar.toAbsolutePath().toString()
        n.normalizeSingularityImageName('foo/baz') == 'docker://foo/baz'
        n.normalizeSingularityImageName('bar.img') == img.toAbsolutePath().toString()

        cleanup:
        base.deleteDir()

    }

    def 'test normalize method for docker' () {
        given:
        def n = Spy(ContainerHandler,constructorArgs:[[engine: 'docker', enabled: true, registry: registry]])

        when:
        def result = n.normalizeImageName(image)

        then:
        1 * n.normalizeDockerImageName(image) >> expected
        result == expected

        where:
        image                       | registry  | expected
        null                        | null      | null
        null                        | 'd.reg'   | null
        'hello'                     | null      | 'hello'
        'cbcrg/hello'               | null      | 'cbcrg/hello'
        'cbcrg/hello'               | 'd.reg'   | 'd.reg/cbcrg/hello'
        'cbcrg/hello'               | 'd.reg/'  | 'd.reg/cbcrg/hello'
        'registry:5000/cbcrg/hello' | 'd.reg'   | 'registry:5000/cbcrg/hello'
    }

    def 'test normalize method for shifter' () {

        given:
        def n = Spy(ContainerHandler,constructorArgs:[[engine: 'shifter', enabled: true]])

        when:
        def result = n.normalizeImageName(image)

        then:
        1 * n.normalizeShifterImageName(image) >> expected
        result == expected

        where:
        image                         | expected
        'busybox'                     | 'docker:busybox:latest'
        'busybox:latest'              | 'docker:busybox:latest'
        'busybox:1.0'                 | 'docker:busybox:1.0'
        'docker:busybox:1.0'          | 'docker:busybox:1.0'
        'hello/world'                 | 'docker:hello/world:latest'
        'hello/world:tag-name'        | 'docker:hello/world:tag-name'
        'docker:hello/world:tag-name' | 'docker:hello/world:tag-name'
        'docker:busybox'              | 'docker:busybox:latest'
    }

    def 'test normalize method for singularity' () {
        given:
        def n = Spy(ContainerHandler,constructorArgs:[[engine: 'singularity', enabled: true]])

        when:
        def result = n.normalizeImageName(image)

        then:
        1 * n.normalizeSingularityImageName(image) >> normalized
        intExpected * n.createCache(n.config, normalized) >> expected
        result == expected

        where:
        image                      | normalized                                       | intExpected | expected
        null                       | null                                             |           0 | null
        ''                         | null                                             |           0 | null
        '/abs/path/bar.img'        | '/abs/path/bar.img'                              |           0 | '/abs/path/bar.img'
        'file:///abs/path/bar.img' | '/abs/path/bar.img'                              |           0 | '/abs/path/bar.img'
        'foo.img'                  | Paths.get('foo.img').toAbsolutePath().toString() |           0 | Paths.get('foo.img').toAbsolutePath().toString()
        'shub://busybox'           | 'shub://busybox'                                 |           1 | '/path/to/busybox'
        'docker://library/busybox' | 'docker://library/busybox'                       |           1 | '/path/to/busybox'
        'foo'                      | 'docker://foo'                                   |           1 | '/path/to/foo'
    }
}
