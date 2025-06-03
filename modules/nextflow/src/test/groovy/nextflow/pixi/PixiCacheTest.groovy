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

package nextflow.pixi

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.util.Duration
import spock.lang.Specification

/**
 *
 * @author Edmund Miller <edmund.miller@seqera.io>
 */
class PixiCacheTest extends Specification {

    def 'should detect TOML file path' () {
        given:
        def cache = new PixiCache()

        expect:
        !cache.isTomlFilePath('python=3.8')
        !cache.isTomlFilePath('env.yaml')
        cache.isTomlFilePath('pixi.toml')
        cache.isTomlFilePath('pyproject.toml')
        cache.isTomlFilePath('env.toml')
        cache.isTomlFilePath('/path/to/pixi.toml')
        !cache.isTomlFilePath('pixi.toml\nsome other content')
    }

    def 'should detect lock file path' () {
        given:
        def cache = new PixiCache()

        expect:
        !cache.isLockFilePath('python=3.8')
        !cache.isLockFilePath('env.yaml')
        !cache.isLockFilePath('pixi.toml')
        cache.isLockFilePath('pixi.lock')
        cache.isLockFilePath('/path/to/pixi.lock')
        !cache.isLockFilePath('pixi.lock\nsome other content')
    }

    def 'should create pixi env prefix path for a package specification' () {
        given:
        def ENV = 'python=3.8'
        def cache = Spy(PixiCache)
        def BASE = Paths.get('/pixi/envs')

        when:
        def prefix = cache.pixiPrefixPath(ENV)
        then:
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/pixi/envs/env-')
        prefix.toString().length() > '/pixi/envs/env-'.length()
    }

    def 'should create pixi env prefix path for a toml file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(PixiCache)
        def BASE = Paths.get('/pixi/envs')
        def ENV = folder.resolve('pixi.toml')
        ENV.text = '''
            [project]
            name = "my-project"
            version = "0.1.0"
            channels = ["conda-forge"]

            [dependencies]
            python = "3.8"
            numpy = ">=1.20"
            '''
            .stripIndent(true)

        when:
        def prefix = cache.pixiPrefixPath(ENV.toString())
        then:
        1 * cache.isTomlFilePath(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/pixi/envs/pixi-')

        cleanup:
        folder?.deleteDir()
    }

    def 'should create pixi env prefix path for a lock file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(PixiCache)
        def BASE = Paths.get('/pixi/envs')
        def ENV = folder.resolve('pixi.lock')
        ENV.text = '''
            version: 3
            environments:
              default:
                channels:
                  - url: https://conda.anaconda.org/conda-forge/
                packages:
                  linux-64:
                    - conda: https://conda.anaconda.org/conda-forge/linux-64/python-3.8.18-h955ad1f_0.conda
            '''
            .stripIndent(true)

        when:
        def prefix = cache.pixiPrefixPath(ENV.toString())
        then:
        1 * cache.isTomlFilePath(ENV.toString())
        1 * cache.isLockFilePath(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/pixi/envs/pixi-')

        cleanup:
        folder?.deleteDir()
    }

    def 'should handle non-existent TOML file' () {
        given:
        def cache = new PixiCache()
        def nonExistentFile = '/non/existent/pixi.toml'

        when:
        cache.pixiPrefixPath(nonExistentFile)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('Pixi environment file does not exist')
    }

    def 'should handle non-existent lock file' () {
        given:
        def cache = new PixiCache()
        def nonExistentFile = '/non/existent/pixi.lock'

        when:
        cache.pixiPrefixPath(nonExistentFile)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('Pixi lock file does not exist')
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
        1 * cache.isLockFilePath(ENV)
        0 * cache.getCacheDir()
        prefix.toString() == folder.toString()

        cleanup:
        folder?.deleteDir()
    }

    def 'should reject prefix path with newlines' () {
        given:
        def cache = new PixiCache()
        def ENV = 'invalid\npath'

        when:
        cache.pixiPrefixPath(ENV)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('Invalid Pixi environment definition')
    }

    def 'should reject non-directory prefix path' () {
        given:
        def cache = new PixiCache()
        def tempFile = Files.createTempFile('test', '.txt')
        def ENV = tempFile.toString()

        when:
        cache.pixiPrefixPath(ENV)

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('Pixi prefix path does not exist or is not a directory')

        cleanup:
        Files.deleteIfExists(tempFile)
    }

    def 'should create a pixi environment for existing prefix' () {
        given:
        def ENV = 'python=3.8'
        def PREFIX = Files.createTempDirectory('foo')
        def cache = Spy(PixiCache)

        when:
        def result = cache.createLocalPixiEnv(ENV, PREFIX)
        then:
        result == PREFIX
        0 * cache.createLocalPixiEnv0(_, _)

        cleanup:
        PREFIX?.deleteDir()
    }

    def 'should create a pixi environment from package specification' () {
        given:
        def ENV = 'python=3.8'
        def PREFIX = Files.createTempDirectory('test-env')
        def cache = Spy(PixiCache)

        when:
        def result = cache.createLocalPixiEnv0(ENV, PREFIX)

        then:
        1 * cache.runCommand(_) >> { String cmd ->
            assert cmd.contains("cd ${PREFIX} && pixi install")
            return 0
        }
        result == PREFIX

        cleanup:
        PREFIX?.deleteDir()
    }

    def 'should create a pixi environment from TOML file' () {
        given:
        def ENV = 'pixi.toml'
        def PREFIX = Files.createTempDirectory('test-env')
        def cache = Spy(PixiCache)

        when:
        def result = cache.createLocalPixiEnv0(ENV, PREFIX)

        then:
        _ * cache.makeAbsolute(ENV) >> Paths.get('/usr/project/pixi.toml')
        1 * cache.runCommand(_) >> { String cmd ->
            assert cmd.contains("pixi install")
            return 0
        }
        result == PREFIX

        cleanup:
        PREFIX?.deleteDir()
    }

    def 'should create pixi env with options' () {
        given:
        def ENV = 'python=3.8'
        def PREFIX = Files.createTempDirectory('test-env')
        def config = new PixiConfig([createOptions: '--verbose'], [:])
        def cache = Spy(new PixiCache(config))

        when:
        def result = cache.createLocalPixiEnv0(ENV, PREFIX)

        then:
        1 * cache.runCommand(_) >> { String cmd ->
            assert cmd.contains("cd ${PREFIX} && pixi install --verbose")
            return 0
        }
        result == PREFIX

        cleanup:
        PREFIX?.deleteDir()
    }

    def 'should get options from the config' () {
        when:
        def cache = new PixiCache(new PixiConfig([:], [:]))
        then:
        cache.createTimeout.minutes == 20
        cache.configCacheDir0 == null
        cache.createOptions == null

        when:
        cache = new PixiCache(new PixiConfig([createTimeout: '5 min', cacheDir: '/pixi/cache', createOptions: '--verbose'], [:]))
        then:
        cache.createTimeout.minutes == 5
        cache.configCacheDir0.toString() == '/pixi/cache'
        cache.createOptions == '--verbose'
    }

    def 'should get cache directory from config' () {
        given:
        def customCacheDir = Files.createTempDirectory('custom-cache')
        def config = new PixiConfig([cacheDir: customCacheDir], [:])
        def cache = Spy(new PixiCache(config))

        when:
        def result = cache.getCacheDir()

        then:
        result.toString() == customCacheDir.toString()

        cleanup:
        customCacheDir?.deleteDir()
    }

    def 'should get cache directory from environment variable' () {
        given:
        def cache = Spy(new PixiCache(new PixiConfig([:], [:])))
        def envCacheDir = Files.createTempDirectory('env-cache')

        when:
        def result = cache.getCacheDir()

        then:
        _ * cache.getEnv() >> [NXF_PIXI_CACHEDIR: envCacheDir.toString()]
        result.toString() == envCacheDir.toString()

        cleanup:
        envCacheDir?.deleteDir()
    }

    def 'should get cache directory from work directory when no config' () {
        given:
        def cache = Spy(new PixiCache(new PixiConfig([:], [:])))
        def workDir = Files.createTempDirectory('work')

        when:
        def result = cache.getCacheDir()

        then:
        1 * cache.getEnv() >> [:]
        1 * cache.getSessionWorkDir() >> workDir
        result.toString() == workDir.resolve('pixi').toString()

        cleanup:
        workDir?.deleteDir()
    }

    def 'should make absolute path' () {
        given:
        def cache = new PixiCache()

        when:
        def result = cache.makeAbsolute('pixi.toml')

        then:
        result.isAbsolute()
        result.toString().endsWith('pixi.toml')
    }

    def 'should handle command execution success' () {
        given:
        def cache = Spy(PixiCache)
        def cmd = 'echo "test"'

        when:
        def result = cache.runCommand(cmd)

        then:
        result == 0
    }

    def 'should handle command execution failure' () {
        given:
        def cache = new PixiCache()
        def cmd = 'false'  // command that always fails

        when:
        cache.runCommand(cmd)

        then:
        def e = thrown(IllegalStateException)
        e.message.contains('Failed to create Pixi environment')
    }

    def 'should handle command execution with timeout' () {
        given:
        def config = new PixiConfig([createTimeout: '1 min'], [:])
        def cache = new PixiCache(config)

        expect:
        cache.createTimeout.minutes == 1
    }

        def 'should get cache path for environment' () {
        given:
        def ENV = 'python=3.8'
        def cache = Spy(PixiCache)
        def BASE = Files.createTempDirectory('pixi-cache')

        when:
        def result = cache.pixiPrefixPath(ENV)

        then:
        1 * cache.getCacheDir() >> BASE
        result.toString().startsWith(BASE.toString())
        result.toString().contains('env-')

        cleanup:
        BASE?.deleteDir()
    }

    def 'should test cache configuration inheritance' () {
        given:
        def CONFIG = [createTimeout: '5 min', cacheDir: '/custom/cache', createOptions: '--verbose']
        def cache = new PixiCache(new PixiConfig(CONFIG, [:]))

        expect:
        cache.createTimeout.minutes == 5
        cache.createOptions == '--verbose'
        cache.configCacheDir0.toString() == '/custom/cache'
    }
}

