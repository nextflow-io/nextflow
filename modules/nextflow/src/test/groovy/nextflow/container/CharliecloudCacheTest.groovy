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

import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Ignore
import spock.lang.Specification
import spock.lang.Timeout
import spock.lang.Unroll
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 * @author Patrick Hüther <patrick.huether@gmail.com>
 */
class CharliecloudCacheTest extends Specification {

    @Unroll
    def 'should return a simple name given an image url'() {

        given:
        def helper = new CharliecloudCache(Mock(ContainerConfig))

        expect:
        helper.simpleName(url) == expected

        where:
        url                         | expected
        'docker://foo/bar'          | 'foo%bar'
        'docker://foo/bar:tag'      | 'foo%bar+tag'
        'shub://hello/world'        | 'hello%world'
        'ftp://hello/world'         | 'hello%world'
        'foo:bar'                   | 'foo+bar'
    }

    def 'should return the cache dir from the config file' () {

        given:
        def dir = Files.createTempDirectory('test')
        and:
        def cacheDir = dir.resolve('nxf.ch')

        when:
        def cache = new CharliecloudCache([cacheDir: "$cacheDir"] as ContainerConfig)
        then:
        cache.getCacheDir() == cacheDir

        cleanup:
        dir.deleteDir()
    }

    def 'should return the cache dir from the environment' () {

        given:
        def dir = Files.createTempDirectory('test')
        and:
        def cacheDir = dir.resolve('nxf.ch')

        when:
        def cache = new CharliecloudCache(GroovyMock(ContainerConfig), [NXF_CHARLIECLOUD_CACHEDIR: "$cacheDir"])
        then:
        cache.getCacheDir() == cacheDir

        cleanup:
        dir.deleteDir()
    }

    def 'should run ch-image pull command'() {

        given:
        def dir = Files.createTempDirectory('test')
        def IMAGE = 'busybox:latest'
        def LOCAL = 'img/busybox+latest'
        def CACHE_PATH = dir.resolve('charliecloud')
        def TARGET_PATH = CACHE_PATH.resolve(LOCAL)
        and:
        def cache = Spy(CharliecloudCache)

        when:
        def result = cache.downloadCharliecloudImage(IMAGE)
        then:
        1 * cache.localImagePath(IMAGE) >> TARGET_PATH
        and:
        1 * cache.runCommand("ch-image pull -s $CACHE_PATH $IMAGE > /dev/null", TARGET_PATH) >> 0
        and:
        CACHE_PATH.parent.exists()
        and:
        result == TARGET_PATH

        cleanup:
        dir.deleteDir()
    }

    def 'should return cached image'() {

        given:
        def dir = Files.createTempDirectory('test')
        def IMAGE = 'busybox:latest'
        def LOCAL = 'img/busybox+latest'
        def CACHE_PATH = dir.resolve('charliecloud')
        def TARGET_PATH = CACHE_PATH.resolve(LOCAL)
        and:
        def cache = Spy(CharliecloudCache)
        TARGET_PATH.mkdirs()

        when:
        def result = cache.downloadCharliecloudImage(IMAGE)
        then:
        1 * cache.localImagePath(IMAGE) >> TARGET_PATH
        0 * cache.runCommand(_) >> 0
        and:
        result == TARGET_PATH

        cleanup:
        dir.deleteDir()
    }

    @Ignore
    @Timeout(1)
    def 'should pull a charliecloud image' () {

        given:
        def IMAGE = 'busybox:latest'
        def LOCAL = 'busybox+latest'
        def dir = Paths.get('/test/path')
        def container = dir.resolve(LOCAL)
        and:
        def cache = Spy(CharliecloudCache)

        when:
        def file = cache.getCachePathFor(IMAGE)
        then:
        1 * cache.downloadCharliecloudImage(IMAGE) >> container
        file == container
    }

}
