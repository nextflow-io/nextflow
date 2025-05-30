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
 *
 */

package nextflow.pixi

import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Specification
/**
 * Tests for PixiCache class
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class PixiCacheTest extends Specification {

    def 'should check if toml file' () {
        given:
        def cache = new PixiCache()

        expect:
        !cache.isTomlFilePath('foo=1.0')
        cache.isTomlFilePath('env.toml')
        cache.isTomlFilePath('pixi.toml')
        cache.isTomlFilePath('foo/bar/env.toml')
    }

    def 'should create pixi env prefix path for a string env' () {
        given:
        def ENV = 'python=3.9'
        def cache = Spy(PixiCache)
        def BASE = Paths.get('/pixi/envs')

        when:
        def prefix = cache.pixiPrefixPath(ENV)
        then:
        1 * cache.isTomlFilePath(ENV)
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == '/pixi/envs/env-7b8a0e5b44da2063a08a3ebfb3ef4be1'
    }

    def 'should create pixi env prefix path for a toml env file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(PixiCache)
        def BASE = Paths.get('/pixi/envs')
        def ENV = folder.resolve('foo.toml')
        def hash = PixiCache.sipHash(ENV)
        ENV.text = '''
            [project]
            name = "foo"
            
            [dependencies]
            python = ">=3.9"
            numpy = "*"
            '''
            .stripIndent(true)  // https://issues.apache.org/jira/browse/GROOVY-9423

        when:
        def prefix = cache.pixiPrefixPath(ENV.toString())
        then:
        1 * cache.isTomlFilePath(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == "/pixi/envs/foo-${hash}"

        cleanup:
        folder?.deleteDir()
    }

    def 'should create pixi env prefix path for a toml env file with name' () {
        given:
        def cache = Spy(PixiCache)
        def BASE = Paths.get('/pixi/envs')
        def ENV = Files.createTempFile('test','.toml')
        def hash = PixiCache.sipHash(ENV)
        ENV.text = '''
            [project]
            name = "my-env-1.1"
            
            [dependencies]
            python = ">=3.9"
            numpy = "*"
            '''
                .stripIndent(true)

        when:
        def prefix = cache.pixiPrefixPath(ENV.toString())
        then:
        1 * cache.isTomlFilePath(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == "/pixi/envs/env-${hash}-myenv11"

        cleanup:
        ENV?.delete()
    }

    def 'should return a pixi prefix directory' () {
        given:
        def cache = Spy(PixiCache)
        def folder = Files.createTempDirectory('test')
        def ENV = folder.toString()

        when:
        def prefix = cache.pixiPrefixPath(ENV)
        then:
        1 * cache.isTomlFilePath(ENV)
        0 * cache.getCacheDir()
        prefix.toString() == folder.toString()

        cleanup:
        folder?.deleteDir()
    }

    def 'should create a pixi environment' () {
        given:
        def ENV = 'python=3.9'
        def PREFIX = Files.createTempDirectory('foo')
        def cache = Spy(PixiCache)

        when:
        // the prefix directory exists ==> no pixi command is executed
        def result = cache.createLocalPixiEnv(ENV, PREFIX)
        then:
        0 * cache.isTomlFilePath(ENV)
        0 * cache.runCommand(_)
        result == PREFIX

        when:
        PREFIX.deleteDir()
        result = cache.createLocalPixiEnv0(ENV, PREFIX)
        then:
        1 * cache.isTomlFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand( "pixi init ${PREFIX} && cd ${PREFIX} && pixi add ${ENV}" ) >> null
        result == PREFIX
    }

    def 'should create a pixi env with a toml file' () {
        given:
        def ENV = 'foo.toml'
        def PREFIX = Paths.get('/pixi/envs/my-env')
        def cache = Spy(PixiCache)

        when:
        def result = cache.createLocalPixiEnv0(ENV, PREFIX)
        then:
        1 * cache.isTomlFilePath(ENV)
        1 * cache.makeAbsolute(ENV) >> Paths.get('/usr/base').resolve(ENV)
        1 * cache.runCommand( "cp /usr/base/foo.toml ${PREFIX}/pixi.toml && pixi install --manifest-path ${PREFIX}/pixi.toml" ) >> null
        result == PREFIX
    }

    def 'should create pixi env with create timeout' () {
        given:
        def ENV = 'python=3.9'
        def PREFIX = Paths.get('/foo/bar')
        and:
        def cache = Spy(new PixiCache(createTimeout: '5 min'))

        when:
        def result = cache.createLocalPixiEnv0(ENV, PREFIX)
        then:
        1 * cache.isTomlFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand( "pixi init ${PREFIX} && cd ${PREFIX} && pixi add ${ENV}" ) >> null
        result == PREFIX
        and:
        cache.createTimeout.minutes == 5
    }

    def 'should get options from the config' () {
        when:
        def cache = new PixiCache(new PixiConfig())
        then:
        cache.createTimeout.minutes == 20
        cache.configCacheDir0 == null

        when:
        cache = new PixiCache(new PixiConfig(createTimeout: '10 min', cacheDir: '/pixi/cache'))
        then:
        cache.createTimeout.minutes == 10
        cache.configCacheDir0 == Paths.get('/pixi/cache')
    }

    def 'should define cache dir from config' () {
        given:
        def folder = Files.createTempDirectory('test'); folder.deleteDir()
        def config = new PixiConfig(cacheDir: folder.toString())
        PixiCache cache = Spy(PixiCache, constructorArgs: [config])

        when:
        def result = cache.getCacheDir()
        then:
        0 * cache.getSessionWorkDir()
        result == folder
        result.exists()

        cleanup:
        folder?.deleteDir()
    }

    def 'should define cache dir from rel path' () {
        given:
        def folder = Paths.get('.test-pixi-cache-' + Math.random())
        def config = new PixiConfig(cacheDir: folder.toString())
        PixiCache cache = Spy(PixiCache, constructorArgs: [config])

        when:
        def result = cache.getCacheDir()
        println result
        then:
        0 * cache.getSessionWorkDir()
        result == folder.toAbsolutePath()
        result.exists()

        cleanup:
        folder?.deleteDir()
    }

    def 'should define cache dir from env' () {
        given:
        def folder = Files.createTempDirectory('test'); folder.deleteDir()
        def config = new PixiConfig()
        PixiCache cache = Spy(PixiCache, constructorArgs: [config])

        when:
        def result = cache.getCacheDir()
        then:
        2 * cache.getEnv() >> [NXF_PIXI_CACHEDIR: folder.toString()]
        0 * cache.getSessionWorkDir()
        result == folder
        result.exists()

        cleanup:
        folder?.deleteDir()
    }

    def 'should define cache dir from session workdir' () {
        given:
        def folder = Files.createTempDirectory('test');
        def cache = Spy(PixiCache)

        when:
        def result = cache.getCacheDir()
        then:
        1 * cache.getSessionWorkDir() >> folder
        result == folder.resolve('pixi')
        result.exists()

        cleanup:
        folder?.deleteDir()
    }
}
