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
import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Ignore
import spock.lang.Specification
import spock.lang.Timeout
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SingularityCacheTest extends Specification {

    @Unroll
    def 'should return a simple name given an image url'() {

        given:
        def helper = new SingularityCache(Mock(ContainerConfig))

        expect:
        helper.simpleName(url) == expected

        where:
        url                         | expected
        'docker://foo/bar'          | 'foo-bar.img'
        'docker://foo/bar:tag'      | 'foo-bar-tag.img'
        'shub://hello/world'        | 'hello-world.img'
        'ftp://hello/world'         | 'hello-world.img'
        'foo:bar'                   | 'foo-bar.img'
        'library://library/default/debian:7'        | 'library-default-debian-7.img'
        'library://library/default/debian.sif'      | 'library-default-debian.sif'
        'library://sylabs/examples/blender.sif:7'   | 'sylabs-examples-blender-7.sif'

    }

    def 'should return the cache dir from the config file' () {

        given:
        def dir = Files.createTempDirectory('test')

        when:
        def cache = new SingularityCache([cacheDir: "$dir"] as ContainerConfig)
        then:
        cache.getCacheDir() == dir

        cleanup:
        dir.deleteDir()
    }

    def 'should return the cache dir from the environment' () {

        given:
        def dir = Files.createTempDirectory('test')

        when:
        def cache = new SingularityCache(Mock(ContainerConfig), [NXF_SINGULARITY_CACHEDIR: "$dir"])
        then:
        cache.getCacheDir() == dir

        cleanup:
        dir.deleteDir()
    }


    def 'should run singularity pull command'() {

        given:
        def dir = Files.createTempDirectory('test')
        def IMAGE = 'docker://pditommaso/foo:latest'
        def LOCAL = 'foo-latest.img'
        ContainerConfig config = [noHttps: true] 
        and:
        def cache = Spy(SingularityCache, constructorArgs: [ config ])

        when:
        cache.downloadSingularityImage(IMAGE)
        then:
        1 * cache.localImagePath(IMAGE) >> dir.resolve(LOCAL)
        1 * cache.runCommand("singularity pull --nohttps --name $LOCAL $IMAGE > /dev/null", dir) >> 0

    }


    def 'should return cached image'() {

        given:
        def dir = Files.createTempDirectory('test')
        def IMAGE = 'docker://pditommaso/foo:latest'
        def LOCAL = 'foo-latest.img'
        def container = dir.resolve(LOCAL)
        container.text = 'dummy'
        and:
        def cache = Spy(SingularityCache)

        when:
        cache.downloadSingularityImage(IMAGE)
        then:
        1 * cache.localImagePath(IMAGE) >> container
        0 * cache.runCommand(_) >> 0

        cleanup:
        dir.deleteDir()
    }

    @Ignore
    @Timeout(1)
    def 'should pull a singularity image' () {

        given:
        def IMAGE = 'docker://pditommaso/foo:latest'
        def LOCAL = 'foo-latest.img'
        def dir = Paths.get('/test/path')
        def container = dir.resolve(LOCAL)
        and:
        def cache = Spy(SingularityCache)

        when:
        def file = cache.getCachePathFor(IMAGE)
        then:
        1 * cache.downloadSingularityImage(IMAGE) >> container
        file == container
    }

}
