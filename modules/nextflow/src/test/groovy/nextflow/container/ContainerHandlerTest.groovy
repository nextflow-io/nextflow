/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
        def n = new ContainerHandler([:], Paths.get('/root/dir'))

        expect:
        n.normalizeSingularityImageName(image) == expected

        where:
        image                      | expected
        null                       | null
        ''                         | null
        '/abs/path/bar.img'        | '/abs/path/bar.img'
        'file:///abs/path/bar.img' | '/abs/path/bar.img'
        'file://foo/bar.img'       | '/root/dir/foo/bar.img'
        'docker://library/busybox' | 'docker://library/busybox'
        'shub://busybox'           | 'shub://busybox'
        'foo://busybox'            | 'foo://busybox'
        'foo'                      | 'docker://foo'
        'foo:2.0'                  | 'docker://foo:2.0'
        'foo.img'                  | 'docker://foo.img'
        'quay.io/busybox'          | 'docker://quay.io/busybox'
        'library://library/default/debian:7'    | 'library://library/default/debian:7'
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

    @Unroll
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

    @Unroll
    def 'test normalize method for singularity' () {
        given:
        def handler = Spy(ContainerHandler,constructorArgs:[[engine: 'singularity', enabled: true]])

        when:
        handler.baseDir = Paths.get('/abs/path/')
        def result = handler.normalizeImageName(IMAGE)

        then:
        1 * handler.normalizeSingularityImageName(IMAGE) >> NORMALIZED
        X * handler.createCache(handler.config, NORMALIZED) >> EXPECTED
        result == EXPECTED

        where:
        IMAGE                                       | NORMALIZED                                        | X           | EXPECTED
        null                                        | null                                              |           0 | null
        ''                                          | null                                              |           0 | null
        '/abs/path/bar.img'                         | '/abs/path/bar.img'                               |           0 | '/abs/path/bar.img'
        '/abs/path bar.img'                         | '/abs/path bar.img'                               |           0 | '/abs/path\\ bar.img'
        'file:///abs/path/bar.img'                  | '/abs/path/bar.img'                               |           0 | '/abs/path/bar.img'
        'foo.img'                                   | Paths.get('foo.img').toAbsolutePath().toString() |       0 | Paths.get('foo.img').toAbsolutePath().toString()
        'shub://busybox'                            | 'shub://busybox'                                  |           1 | '/path/to/busybox'
        'docker://library/busybox'                  | 'docker://library/busybox'                        |           1 | '/path/to/busybox'
        'foo'                                       | 'docker://foo'                                    |           1 | '/path/to/foo'
        'library://pditommaso/foo/bar.sif:latest'   | 'library://pditommaso/foo/bar.sif:latest'         |           1 | '/path/to/foo-bar-latest.img'
    }

    def 'should not invoke caching when engine is disabled' () {
        given:
        final handler = Spy(ContainerHandler,constructorArgs:[[engine: 'singularity']])
        final IMAGE = 'docker://foo.img'

        when:
        handler.config.enabled = false
        def result = handler.normalizeImageName(IMAGE)
        then:
        1 * handler.normalizeSingularityImageName(IMAGE) >> IMAGE
        0 * handler.createCache(_,_) >> null
        result == IMAGE

        when:
        handler.config.enabled = true
        result = handler.normalizeImageName(IMAGE)
        then:
        1 * handler.normalizeSingularityImageName(IMAGE) >> IMAGE
        1 * handler.createCache(_,IMAGE) >> '/some/path/foo.img'
        result == '/some/path/foo.img'
    }
}
