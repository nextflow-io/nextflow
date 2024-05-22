/*
 * Copyright 2013-2024, Seqera Labs
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package nextflow.container

import nextflow.executor.Executor
import spock.lang.Specification

import java.nio.file.Paths
import java.nio.file.Files

import spock.lang.Unroll

/**
 * @author Emilio Palumbo <emilio.palumbo@crg.eu>
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

    @Unroll
    def 'test normalize singularity image #image' () {

        given:
        def n = new ContainerHandler([registry: registry], Paths.get('/root/dir'))

        expect:
        n.normalizeSingularityImageName(image) == expected

        where:
        image                      | registry   | expected
        null                       | null       | null
        ''                         | null       | null
        '/abs/path/bar.img'        | null       | '/abs/path/bar.img'
        'file:///abs/path/bar.img' | null       | '/abs/path/bar.img'
        'file://foo/bar.img'       | null       | '/root/dir/foo/bar.img'
        'docker://library/busybox' | null       | 'docker://library/busybox'
        'shub://busybox'           | null       | 'shub://busybox'
        'foo://busybox'            | null       | 'foo://busybox'
        'foo'                      | null       | 'docker://foo'
        'foo:2.0'                  | null       | 'docker://foo:2.0'
        'foo.img'                  | null       | 'docker://foo.img'
        'quay.io/busybox'          | null       | 'docker://quay.io/busybox'
        'library://library/default/debian:7'    | null       | 'library://library/default/debian:7'
        'http://reg.io/v1/alpine:latest'        | null       | 'http://reg.io/v1/alpine:latest'
        'https://reg.io/v1/alpine:latest'       | null       | 'https://reg.io/v1/alpine:latest'
        and:
        '/abs/path/bar.img'        | 'my.reg'  | '/abs/path/bar.img'
        'quay.io/busybox'          | 'my.reg'  | 'docker://quay.io/busybox'
        'foo'                      | 'my.reg'  | 'docker://my.reg/foo'
        'foo:2.0'                  | 'my.reg'  | 'docker://my.reg/foo:2.0'
        'foo.img'                  | 'my.reg'  | 'docker://my.reg/foo.img'
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

    @Unroll
    def 'test normalize method for docker' () {
        given:
        def n = Spy(new ContainerHandler([engine: 'docker', enabled: true, registry: registry]))

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

    @Unroll
    def 'test normalize method for shifter' () {

        given:
        def n = Spy(new ContainerHandler([engine: 'shifter', enabled: true]))

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

    def 'should use docker for container native'  () {
        given:
        def EXECUTOR  = Mock(Executor)
        def IMAGE = 'foo:latest'
        def handler = Spy(new ContainerHandler([engine: 'shifter', enabled: true]))

        when:
        def result = handler.normalizeImageName(IMAGE)
        then:
        1 * handler.normalizeShifterImageName(IMAGE) >> 'shifter://image'
        0 * handler.normalizeDockerImageName(IMAGE) >> null
        and:
        result == 'shifter://image'
            
    }

    @Unroll
    def 'test normalize method for charliecloud' () {

       given:
        def n = new ContainerHandler([registry: registry])

        expect:
        n.normalizeCharliecloudImageName(image) == expected

        where:
        image                      | registry   | expected
        null                       | null       | null
        ''                         | null       | null
        '/abs/path/bar.img'        | null       | '/abs/path/bar.img'
        'docker://library/busybox' | null       | 'library/busybox:latest'
        'shub://busybox'           | null       | 'shub://busybox'
        'foo://busybox'            | null       | 'foo://busybox'
        'foo'                      | null       | 'foo:latest'
        'foo:2.0'                  | null       | 'foo:2.0'
        'foo.img'                  | null       | 'foo.img:latest'
        'quay.io/busybox'          | null       | 'quay.io/busybox:latest'
        'http://reg.io/v1/alpine:latest'        | null       | 'http://reg.io/v1/alpine:latest'
        'https://reg.io/v1/alpine:latest'       | null       | 'https://reg.io/v1/alpine:latest'
        and:
        '/abs/path/bar.img'        | 'my.reg'  | '/abs/path/bar.img'
        'busybox'                  | 'my.reg'  | 'my.reg/busybox:latest'
        'foo:2.0'                  | 'my.reg'  | 'my.reg/foo:2.0'
    }

    @Unroll
    def 'test normalize method for singularity' () {
        given:
        def BASE = Paths.get('/abs/path/')
        def handler = Spy(new ContainerHandler(engine: 'singularity', enabled: true, ociMode: OCI, baseDir: BASE))

        when:
        def result = handler.normalizeImageName(IMAGE)

        then:
        1 * handler.normalizeSingularityImageName(IMAGE) >> NORMALIZED
        X * handler.createSingularityCache(handler.config, NORMALIZED) >> EXPECTED
        result == EXPECTED

        where:
        IMAGE                                       | NORMALIZED                                        | OCI   | X           | EXPECTED
        null                                        | null                                              | false |           0 | null
        ''                                          | null                                              | false |           0 | null
        '/abs/path/bar.img'                         | '/abs/path/bar.img'                               | false |           0 | '/abs/path/bar.img'
        '/abs/path bar.img'                         | '/abs/path bar.img'                               | false |           0 | '/abs/path\\ bar.img'
        'file:///abs/path/bar.img'                  | '/abs/path/bar.img'                               | false |           0 | '/abs/path/bar.img'
        'foo.img'                                   | Paths.get('foo.img').toAbsolutePath().toString()  | false |           0 | Paths.get('foo.img').toAbsolutePath().toString()
        'shub://busybox'                            | 'shub://busybox'                                  | false |           1 | '/path/to/busybox'
        'docker://library/busybox'                  | 'docker://library/busybox'                        | false |           1 | '/path/to/busybox'
        'foo'                                       | 'docker://foo'                                    | false |           1 | '/path/to/foo'
        'library://pditommaso/foo/bar.sif:latest'   | 'library://pditommaso/foo/bar.sif:latest'         | false |           1 | '/path/to/foo-bar-latest.img'
        and:
        'docker://library/busybox'                  | 'docker://library/busybox'                        | true  |           0 | 'docker://library/busybox'
        'shub://busybox'                            | 'shub://busybox'                                  | true  |           1 | '/path/to/busybox'

    }

    @Unroll
    def 'test normalize method for OCI direct mode' () {
        given:
        def BASE = Paths.get('/abs/path/')
        def handler = Spy(new ContainerHandler(engine: 'apptainer', enabled: true, ociAutoPull:AUTO, baseDir: BASE))

        when:
        def result = handler.normalizeImageName(IMAGE)

        then:
        1 * handler.normalizeApptainerImageName(IMAGE) >> NORMALIZED
        X * handler.createApptainerCache(handler.config, NORMALIZED) >> EXPECTED
        result == EXPECTED

        where:
        IMAGE                                       | NORMALIZED                                        | ENGINE            | AUTO  | X | EXPECTED
        null                                        | null                                              | 'singularity'     | false |           0 | null
        ''                                          | null                                              | 'singularity'     | false |           0 | null
        '/abs/path/bar.img'                         | '/abs/path/bar.img'                               | 'singularity'     | false |           0 | '/abs/path/bar.img'
        '/abs/path bar.img'                         | '/abs/path bar.img'                               | 'singularity'     | false |           0 | '/abs/path\\ bar.img'
        'file:///abs/path/bar.img'                  | '/abs/path/bar.img'                               | 'singularity'     | false |           0 | '/abs/path/bar.img'
        'foo.img'                                   | Paths.get('foo.img').toAbsolutePath().toString()  | 'singularity'     | false |           0 | Paths.get('foo.img').toAbsolutePath().toString()
        'shub://busybox'                            | 'shub://busybox'                                  | 'singularity'     | false |           1 | '/path/to/busybox'
        'docker://library/busybox'                  | 'docker://library/busybox'                        | 'singularity'     | false |           1 | '/path/to/busybox'
        'foo'                                       | 'docker://foo'                                    | 'singularity'     | false |           1 | '/path/to/foo'
        'library://pditommaso/foo/bar.sif:latest'   | 'library://pditommaso/foo/bar.sif:latest'         | 'singularity'     | false |           1 | '/path/to/foo-bar-latest.img'
        and:
        'docker://library/busybox'                  | 'docker://library/busybox'                        | 'singularity'     | true  |           0 | 'docker://library/busybox'
        'shub://busybox'                            | 'shub://busybox'                                  | 'singularity'     | true  |           1 | '/path/to/busybox'

        and:
        null                                        | null                                              | 'apptainer'     | false |           0 | null
        ''                                          | null                                              | 'apptainer'     | false |           0 | null
        '/abs/path/bar.img'                         | '/abs/path/bar.img'                               | 'apptainer'     | false |           0 | '/abs/path/bar.img'
        '/abs/path bar.img'                         | '/abs/path bar.img'                               | 'apptainer'     | false |           0 | '/abs/path\\ bar.img'
        'file:///abs/path/bar.img'                  | '/abs/path/bar.img'                               | 'apptainer'     | false |           0 | '/abs/path/bar.img'
        'foo.img'                                   | Paths.get('foo.img').toAbsolutePath().toString()  | 'apptainer'     | false |           0 | Paths.get('foo.img').toAbsolutePath().toString()
        'shub://busybox'                            | 'shub://busybox'                                  | 'apptainer'     | false |           1 | '/path/to/busybox'
        'docker://library/busybox'                  | 'docker://library/busybox'                        | 'apptainer'     | false |           1 | '/path/to/busybox'
        'foo'                                       | 'docker://foo'                                    | 'apptainer'     | false |           1 | '/path/to/foo'
        'library://pditommaso/foo/bar.sif:latest'   | 'library://pditommaso/foo/bar.sif:latest'         | 'apptainer'     | false |           1 | '/path/to/foo-bar-latest.img'
        and:
        'docker://library/busybox'                  | 'docker://library/busybox'                        | 'apptainer'     | true  |           0 | 'docker://library/busybox'
        'shub://busybox'                            | 'shub://busybox'                                  | 'apptainer'     | true  |           1 | '/path/to/busybox'


    }

    def 'should not invoke caching when engine is disabled' () {
        given:
        final handler = Spy(new ContainerHandler([engine: 'singularity']))
        final IMAGE = 'docker://foo.img'

        when:
        handler.config.enabled = false
        def result = handler.normalizeImageName(IMAGE)
        then:
        1 * handler.normalizeSingularityImageName(IMAGE) >> IMAGE
        0 * handler.createSingularityCache(_,_) >> null
        result == IMAGE

        when:
        handler.config.enabled = true
        result = handler.normalizeImageName(IMAGE)
        then:
        1 * handler.normalizeSingularityImageName(IMAGE) >> IMAGE
        1 * handler.createSingularityCache(_,IMAGE) >> '/some/path/foo.img'
        result == '/some/path/foo.img'
    }

    def 'should invoke singularity cache' () {
        given:
        def handler = Spy(ContainerHandler,constructorArgs:[[engine: 'singularity', enabled: true]])

        when:
        def result = handler.normalizeImageName(IMG)
        then:
        TIMES * handler.createSingularityCache(_, NORM) >> EXPECTED
        
        then:
        result == EXPECTED

        where:
        IMG                     | NORM                  | TIMES | EXPECTED
        'foo'                   | 'docker://foo'        | 1     | '/local/img/foo'
        'library://foo:latest'  | 'library://foo:latest'| 1     | '/local/img/foo.img'
        'http://bar:latest'     | 'http://bar:latest'   | 1     | '/local/http/foo.img'
        'https://bar:latest'    | 'https://bar:latest'  | 1     | '/local/https/foo.img'
        '/some/container.img'   | '/some/container.img' | 0     | '/some/container.img'
    }

}
