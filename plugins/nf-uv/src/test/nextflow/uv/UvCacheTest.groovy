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

package nextflow.uv

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.SysEnv
import spock.lang.Specification

class UvCacheTest extends Specification {

    def setupSpec() {
        SysEnv.push([:])
    }

    def cleanupSpec() {
        SysEnv.pop()
    }

    def 'should detect requirements file' () {
        given:
        def cache = new UvCache()

        expect:
        !cache.isRequirementsFile('numpy')
        cache.isRequirementsFile('requirements.txt')
        cache.isRequirementsFile('requirements.in')
        !cache.isRequirementsFile("requirements.txt\nfoo")
    }

    def 'should detect pyproject file' () {
        given:
        def cache = new UvCache()

        expect:
        !cache.isPyProjectFile('numpy')
        cache.isPyProjectFile('pyproject.toml')
        cache.isPyProjectFile('/path/to/pyproject.toml')
        !cache.isPyProjectFile("pyproject.toml\nfoo")
    }

    def 'should create uv env prefix path for a string env' () {
        given:
        def ENV = 'numpy pandas'
        def cache = Spy(UvCache)
        def BASE = Paths.get('/uv/envs')

        when:
        def prefix = cache.uvPrefixPath(ENV)
        then:
        1 * cache.isRequirementsFile(ENV)
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/uv/envs/env-')
    }

    def 'should create uv env prefix path for a requirements file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(UvCache)
        def BASE = Paths.get('/uv/envs')
        def ENV = folder.resolve('requirements.txt')
        ENV.text = '''\
            numpy==1.24.0
            pandas>=2.0
            '''.stripIndent()

        when:
        def prefix = cache.uvPrefixPath(ENV.toString())
        then:
        1 * cache.isRequirementsFile(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/uv/envs/env-')

        cleanup:
        folder?.deleteDir()
    }

    def 'should return existing directory for path with slash' () {
        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(UvCache)

        when:
        def prefix = cache.uvPrefixPath(folder.toString())
        then:
        prefix == folder

        cleanup:
        folder?.deleteDir()
    }

    def 'should throw for non-existing directory path' () {
        given:
        def cache = Spy(UvCache)
        def ENV = '/non/existing/path'

        when:
        cache.uvPrefixPath(ENV)
        then:
        thrown(IllegalArgumentException)
    }

    def 'should throw for multi-line env' () {
        given:
        def cache = Spy(UvCache)
        def ENV = "numpy\npandas"

        when:
        cache.uvPrefixPath(ENV)
        then:
        thrown(IllegalArgumentException)
    }

    def 'should create the correct uv venv command for package list' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-abc123')
        def cache = Spy(UvCache)
        cache.@installOptions = null
        cache.@pythonVersion = null
        cache.@createTimeout = nextflow.util.Duration.of('20min')

        when:
        cache.createLocalUvEnv0('numpy pandas', prefixPath)
        then:
        1 * cache.runCommand({ String cmd ->
            cmd.contains('uv venv') && cmd.contains('uv pip install') && cmd.contains('numpy pandas')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should create the correct uv venv command with python version' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-abc123')
        def cache = Spy(UvCache)
        cache.@installOptions = null
        cache.@pythonVersion = '3.11'
        cache.@createTimeout = nextflow.util.Duration.of('20min')

        when:
        cache.createLocalUvEnv0('numpy', prefixPath)
        then:
        1 * cache.runCommand({ String cmd ->
            cmd.contains('uv venv --python 3.11') && cmd.contains('uv pip install')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should create the correct uv venv command for requirements file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def reqFile = folder.resolve('requirements.txt')
        reqFile.text = 'numpy==1.24.0\npandas>=2.0'
        def prefixPath = folder.resolve('env-abc123')
        def cache = Spy(UvCache)
        cache.@installOptions = null
        cache.@pythonVersion = null
        cache.@createTimeout = nextflow.util.Duration.of('20min')

        when:
        cache.createLocalUvEnv0(reqFile.toString(), prefixPath)
        then:
        1 * cache.runCommand({ String cmd ->
            cmd.contains('uv venv') && cmd.contains('-r') && cmd.contains('requirements.txt')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should create the correct uv venv command for a custom-named requirements file' () {
        given:
        def folder = Files.createTempDirectory('test')
        // non-conventional name but supported `.in` extension
        def reqFile = folder.resolve('deps.in')
        reqFile.text = 'numpy==1.24.0\npandas>=2.0'
        def prefixPath = folder.resolve('env-abc123')
        def cache = Spy(UvCache)
        cache.@installOptions = null
        cache.@pythonVersion = null
        cache.@createTimeout = nextflow.util.Duration.of('20min')

        when:
        cache.createLocalUvEnv0(reqFile.toString(), prefixPath)
        then:
        1 * cache.runCommand({ String cmd ->
            cmd.contains('uv pip install') && cmd.contains('-r') && cmd.contains('deps.in')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should create the correct uv venv command for a pyproject file' () {
        given:
        def folder = Files.createTempDirectory('test')
        def pyprojectFile = folder.resolve('pyproject.toml')
        pyprojectFile.text = '[project]\nname = "demo"\ndependencies = ["numpy"]\n'
        def prefixPath = folder.resolve('env-abc123')
        def cache = Spy(UvCache)
        cache.@installOptions = null
        cache.@pythonVersion = null
        cache.@createTimeout = nextflow.util.Duration.of('20min')

        when:
        cache.createLocalUvEnv0(pyprojectFile.toString(), prefixPath)
        then:
        1 * cache.runCommand({ String cmd ->
            cmd.contains('uv venv') && cmd.contains('uv pip install') && cmd.contains('-r') && cmd.contains('pyproject.toml')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should include python version in hash' () {
        given:
        def cache1 = Spy(UvCache)
        cache1.@pythonVersion = null
        def cache2 = Spy(UvCache)
        cache2.@pythonVersion = '3.12'
        def BASE = Paths.get('/uv/envs')

        when:
        def prefix1 = cache1.uvPrefixPath('numpy')
        def prefix2 = cache2.uvPrefixPath('numpy')
        then:
        _ * cache1.getCacheDir() >> BASE
        _ * cache2.getCacheDir() >> BASE
        prefix1 != prefix2
    }

    def 'should apply a per-process install-options override in the command' () {
        given:
        def folder = Files.createTempDirectory('test')
        def prefixPath = folder.resolve('env-ovr')
        def cache = Spy(UvCache)
        cache.@installOptions = '--from-config'   // config-level default
        cache.@pythonVersion = null
        cache.@createTimeout = nextflow.util.Duration.of('20min')

        when:
        cache.createLocalUvEnv0('numpy', prefixPath, '--no-cache')
        then:
        1 * cache.runCommand({ String cmd ->
            cmd.contains('--no-cache') && !cmd.contains('--from-config')
        }) >> 0

        cleanup:
        folder?.deleteDir()
    }

    def 'should include the install-options override in the env hash' () {
        given:
        def base = Spy(UvCache); base.@installOptions = null; base.@pythonVersion = null
        def ovr  = Spy(UvCache); ovr.@installOptions = null;  ovr.@pythonVersion = null
        def BASE = Paths.get('/uv/envs')

        when:
        def p1 = base.uvPrefixPath('numpy')
        def p2 = ovr.uvPrefixPath('numpy', '--no-cache')
        then:
        _ * base.getCacheDir() >> BASE
        _ * ovr.getCacheDir() >> BASE
        p1 != p2
    }

}
