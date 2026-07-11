/*
 * Copyright 2013-2026, Seqera Labs
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

import spock.lang.Specification

class PixiCacheTest extends Specification {

    def 'should apply a per-process create-options override in the command' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-ovr')
        def cache = Spy(PixiCache)
        cache.@createOptions = '--from-config'   // config-level default
        cache.@createTimeout = nextflow.util.Duration.of('20min')

        when:
        cache.createLocalPixiEnv0('cowpy', prefixPath, '--frozen')
        then:
        1 * cache.runCommand({ String cmd ->
            cmd.contains('pixi install') && cmd.contains('--frozen') && !cmd.contains('--from-config')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should include the create-options override in the env hash' () {
        given:
        def base = Spy(PixiCache); base.@createOptions = null
        def ovr  = Spy(PixiCache); ovr.@createOptions = null
        def BASE = Paths.get('/pixi/envs')

        when:
        def p1 = base.pixiPrefixPath('cowpy')
        def p2 = ovr.pixiPrefixPath('cowpy', '--frozen')
        then:
        _ * base.getCacheDir() >> BASE
        _ * ovr.getCacheDir() >> BASE
        p1 != p2
    }

    def 'should build an install command for a multi-package spec' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-multi')
        def cache = Spy(PixiCache)
        cache.@createOptions = null
        cache.@createTimeout = nextflow.util.Duration.of('20min')

        when:
        cache.createLocalPixiEnv0('cowpy numpy>=1.20', prefixPath)
        then:
        1 * cache.runCommand({ String cmd -> cmd.contains('pixi install') }) >> 0
        and:
        // a generated pixi.toml with the requested deps is written into the prefix
        def manifest = prefixPath.resolve('pixi.toml').text
        manifest.contains('cowpy = "*"')
        manifest.contains('numpy = ">=1.20"')

        cleanup:
        folder?.deleteDir()
    }

    def 'should resolve a prefix for a custom-named toml file and read its content' () {
        given:
        def folder = Files.createTempDirectory('test')
        def tomlFile = folder.resolve('myenv.toml')
        tomlFile.text = '''\
            [dependencies]
            cowpy = "*"
            '''.stripIndent()
        def cache = Spy(PixiCache)
        cache.@createOptions = null
        def BASE = folder.resolve('cache')

        when:
        def prefix = cache.pixiPrefixPath(tomlFile.toString())
        then:
        _ * cache.getCacheDir() >> BASE
        and:
        // baseName of the custom file drives the prefix name
        prefix.name.startsWith('myenv-')
        prefix.parent == BASE

        cleanup:
        folder?.deleteDir()
    }

    def 'should build a pixi install command for a custom-named toml file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def tomlFile = folder.resolve('myenv.toml')
        tomlFile.text = '''\
            [dependencies]
            cowpy = "*"
            '''.stripIndent()
        def prefixPath = folder.resolve('env-toml')
        def cache = Spy(PixiCache)
        cache.@createOptions = null
        cache.@createTimeout = nextflow.util.Duration.of('20min')

        when:
        cache.createLocalPixiEnv0(tomlFile.toString(), prefixPath)
        then:
        1 * cache.runCommand({ String cmd ->
            cmd.contains('pixi install') && cmd.contains(folder.toString())
        }) >> 0
        and:
        // the .pixi marker file points back to the project dir
        def marker = prefixPath.resolve('.pixi')
        Files.isRegularFile(marker)
        marker.text == folder.toString()

        cleanup:
        folder?.deleteDir()
    }

    def 'should throw when the toml file does not exist' () {
        given:
        def folder = Files.createTempDirectory('test')
        def missing = folder.resolve('does-not-exist.toml')
        def cache = Spy(PixiCache)
        cache.@createOptions = null

        when:
        cache.pixiPrefixPath(missing.toString())
        then:
        thrown(IllegalArgumentException)

        cleanup:
        folder?.deleteDir()
    }

    def 'should throw when a slash path is not an existing directory' () {
        given:
        def folder = Files.createTempDirectory('test')
        def missingDir = folder.resolve('no/such/dir')
        def cache = Spy(PixiCache)
        cache.@createOptions = null

        when:
        cache.pixiPrefixPath(missingDir.toString())
        then:
        thrown(IllegalArgumentException)

        cleanup:
        folder?.deleteDir()
    }

    def 'should throw for a multi-line/invalid environment definition' () {
        given:
        def cache = Spy(PixiCache)
        cache.@createOptions = null

        when:
        cache.pixiPrefixPath('cowpy\nnumpy')
        then:
        thrown(IllegalArgumentException)
    }
}
