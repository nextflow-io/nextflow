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

package nextflow.conda

import java.nio.file.Files
import java.nio.file.Paths

import nextflow.SysEnv
import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CondaCacheTest extends Specification {

    def setupSpec() {
        SysEnv.push([:])
    }

    def cleanupSpec() {
        SysEnv.pop()
    }

    def 'should env file' () {

        given:
        def cache = new CondaCache()

        expect:
        !cache.isYamlFilePath('foo=1.0')
        cache.isYamlFilePath('env.yml')
        cache.isYamlFilePath('env.yaml')
    }

    def 'should text file' () {

        given:
        def cache = new CondaCache()

        expect:
        !cache.isTextFilePath('foo=1.0')
        !cache.isTextFilePath('env.yaml')
        !cache.isTextFilePath('foo.txt\nbar.txt')
        cache.isTextFilePath('env.txt')
        cache.isTextFilePath('foo/bar/env.txt')
    }

    def 'should detect lock file content' () {

        given:
        def cache = new CondaCache()

        expect:
        // Valid lock file content - @EXPLICIT marker present
        cache.isLockFile('@EXPLICIT\nhttps://conda.anaconda.org/conda-forge/linux-64/package-1.0.tar.bz2')
        cache.isLockFile('# This file may be used to create an environment\n@EXPLICIT\nhttps://url')
        cache.isLockFile('# comment\n# another comment\n@EXPLICIT\nhttps://url')
        // With spaces/indentation
        cache.isLockFile('  @EXPLICIT  \nhttps://url')

        // Invalid lock file content - no @EXPLICIT marker
        !cache.isLockFile('foo=1.0')
        !cache.isLockFile('channels:\n  - conda-forge\ndependencies:\n  - bwa')
        !cache.isLockFile('')
        !cache.isLockFile(null)
        // @EXPLICIT after 20 lines should not be detected
        !cache.isLockFile((1..25).collect { "# line $it" }.join('\n') + '\n@EXPLICIT')
    }

    def 'should detect remote file' () {

        given:
        def cache = new CondaCache()

        expect:
        // HTTP/HTTPS protocols
        cache.isRemoteFile('http://example.com/condalock')
        cache.isRemoteFile('https://example.com/env.lock')
        cache.isRemoteFile('https://wave.seqera.io/v1alpha1/builds/bd-123/condalock')
        // Cloud storage protocols
        cache.isRemoteFile('s3://bucket/path/to/condalock')
        cache.isRemoteFile('gs://bucket/path/to/condalock')
        cache.isRemoteFile('az://container/path/to/condalock')
        // FTP protocol
        cache.isRemoteFile('ftp://example.com/file')
        // Not remote files
        !cache.isRemoteFile('foo.yml')
        !cache.isRemoteFile('/path/to/env.lock')
        !cache.isRemoteFile('file:///path/to/env.lock')
        !cache.isRemoteFile(null)
    }

    def 'should create conda env prefix path for a string env' () {

        given:
        def ENV = 'bwa=1.7.2'
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')

        when:
        def prefix = cache.condaPrefixPath(ENV)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == '/conda/envs/env-eaeb133f4ca62c95e9c0eec7ef8d553b'
    }

    def 'should create conda env prefix path for remote uri' () {

        given:
        def ENV = 'https://foo.com/lock-file.yml'
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')
        def STAGED_CONTENT = 'name: test\ndependencies:\n  - bwa'

        when:
        def prefix = cache.condaPrefixPath(ENV)
        then:
        1 * cache.isRemoteFile(ENV) >> true
        1 * cache.stageRemoteFile(ENV) >> {
            def tempFile = Files.createTempFile('test', '.yml')
            tempFile.text = STAGED_CONTENT
            return tempFile
        }
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith('/conda/envs/env-')
    }

    def 'should create conda env prefix path for a yaml env file' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')
        def ENV = folder.resolve('foo.yml')
        ENV.text = '''
            channels:
              - conda-forge
              - bioconda
            dependencies:
              # Default bismark
              - star=2.5.4a
              - bwa=0.7.15
            '''
            .stripIndent(true)  // https://issues.apache.org/jira/browse/GROOVY-9423
        when:
        def prefix = cache.condaPrefixPath(ENV.toString())
        then:
        1 * cache.isYamlFilePath(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == "/conda/envs/env-64874f9dc9e7be788384bccef357a4f4"

        cleanup:
        folder?.deleteDir()

    }

    def 'should create conda env prefix path for a env yaml file with name' () {

        given:
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')
        def ENV = Files.createTempFile('test','.yml')
        ENV.text = '''
            name: my-env-1.1
            channels:
              - conda-forge
              - bioconda
            dependencies:
              # Default bismark
              - star=2.5.4a
              - bwa=0.7.15
            '''
                .stripIndent(true)

        when:
        def prefix = cache.condaPrefixPath(ENV.toString())
        then:
        1 * cache.isYamlFilePath(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == "/conda/envs/env-5b5c72e839d0c7dcabb5d06607c205fc"

    }

    def 'should create conda env prefix path for a text env file' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')
        def ENV = folder.resolve('bar.txt')
        ENV.text = '''
                star=2.5.4a
                bwa=0.7.15
                multiqc=1.2.3
                '''
                .stripIndent(true)  // https://issues.apache.org/jira/browse/GROOVY-9423

        when:
        def prefix = cache.condaPrefixPath(ENV.toString())
        then:
        1 * cache.isYamlFilePath(ENV.toString())
        1 * cache.isTextFilePath(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == "/conda/envs/env-85371202d8820331ff19ae89c0595497"

        cleanup:
        folder?.deleteDir()

    }

    def 'should return a conda prefix directory' () {

        given:
        def cache = Spy(CondaCache)
        def folder = Files.createTempDirectory('test')
        def ENV = folder.toString()

        when:
        def prefix = cache.condaPrefixPath(ENV)
        then:
        1 * cache.isYamlFilePath(ENV)
        0 * cache.getCacheDir()
        prefix.toString() == folder.toString()

        cleanup:
        folder?.deleteDir()
    }

    def 'should create conda env prefix path for a lock file with .lock extension' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')
        def ENV = folder.resolve('env.lock')
        ENV.text = '''
                # This file may be used to create an environment using:
                # $ conda create --name <env> --file <this file>
                @EXPLICIT
                https://conda.anaconda.org/conda-forge/linux-64/package-1.0.tar.bz2
                '''
                .stripIndent(true)

        when:
        def prefix = cache.condaPrefixPath(ENV.toString())
        then:
        1 * cache.isYamlFilePath(ENV.toString())
        1 * cache.isTextFilePath(ENV.toString())
        1 * cache.isLockFile(_) >> true
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith("/conda/envs/env-")

        cleanup:
        folder?.deleteDir()
    }

    def 'should create conda env prefix path for a lock file with no extension' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')
        def ENV = folder.resolve('condalock')
        ENV.text = '''
                @EXPLICIT
                https://conda.anaconda.org/conda-forge/linux-64/package-1.0.tar.bz2
                '''
                .stripIndent(true)

        when:
        def prefix = cache.condaPrefixPath(ENV.toString())
        then:
        1 * cache.isYamlFilePath(ENV.toString())
        1 * cache.isTextFilePath(ENV.toString())
        1 * cache.isLockFile(_) >> true
        1 * cache.getCacheDir() >> BASE
        prefix.toString().startsWith("/conda/envs/env-")

        cleanup:
        folder?.deleteDir()
    }

    def 'should reject file with non-standard extension that is not a lock file' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(CondaCache)
        def ENV = folder.resolve('env.lock')
        ENV.text = '''
                # This is not a valid lock file
                bwa=1.0
                samtools=1.15
                '''
                .stripIndent(true)

        when:
        cache.condaPrefixPath(ENV.toString())

        then:
        def e = thrown(IllegalArgumentException)
        e.message.contains('not a valid lock file')

        cleanup:
        folder?.deleteDir()
    }

    def 'should create a conda environment' () {
        given:
        def ENV = 'bwa=1.1.1'
        def PREFIX = Files.createTempDirectory('foo')
        def cache = Spy(CondaCache)

        when:
        // the prefix directory exists ==> no conda command is executed
        def result = cache.createLocalCondaEnv(ENV, PREFIX)
        then:
        0 * cache.isYamlFilePath(ENV)
        0 * cache.runCommand(_)
        result == PREFIX

        when:
        PREFIX.deleteDir()
        result = cache.createLocalCondaEnv0(ENV,PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand( "conda create --yes --quiet --prefix $PREFIX $ENV" ) >> null
        result == PREFIX
    }

    def 'should create a conda environment - using mamba' () {
        given:
        def ENV = 'bwa=1.1.1'
        def PREFIX = Files.createTempDirectory('foo')
        def cache = Spy(new CondaCache(useMamba: true))

        when:
        // the prefix directory exists ==> no mamba command is executed
        def result = cache.createLocalCondaEnv(ENV, PREFIX)
        then:
        0 * cache.isYamlFilePath(ENV)
        0 * cache.runCommand(_)
        result == PREFIX

        when:
        PREFIX.deleteDir()
        result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand("mamba create --yes --quiet --prefix $PREFIX $ENV") >> null
        result == PREFIX
    }

    def 'should create a conda environment - using micromamba' () {
        given:
        def ENV = 'bwa=1.1.1'
        def PREFIX = Files.createTempDirectory('foo')
        def cache = Spy(new CondaCache(useMicromamba: true))

        when:
        // the prefix directory exists ==> no mamba command is executed
        def result = cache.createLocalCondaEnv(ENV, PREFIX)
        then:
        0 * cache.isYamlFilePath(ENV)
        0 * cache.runCommand(_)
        result == PREFIX

        when:
        PREFIX.deleteDir()
        result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand("micromamba create --yes --quiet --prefix $PREFIX $ENV") >> null
        result == PREFIX
    }

    def 'should create a conda environment using mamba and remote yaml file' () {
        given:
        def ENV = 'http://foo.com/some/env.yml'
        def PREFIX = Files.createTempDirectory('foo')
        def STAGED_PATH = Paths.get('/staged/env.yml')
        def cache = Spy(new CondaCache(useMamba: true))

        when:
        // the prefix directory exists ==> no mamba command is executed
        def result = cache.createLocalCondaEnv(ENV, PREFIX)
        then:
        0 * cache.isYamlFilePath(ENV)
        0 * cache.runCommand(_)
        result == PREFIX

        when:
        PREFIX.deleteDir()
        result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        // isYamlFilePath returns true for .yml files, so we enter YAML branch
        1 * cache.isYamlFilePath(ENV)
        // getLocalFilePath is called, which internally calls isRemoteFile and stageRemoteFile
        1 * cache.getLocalFilePath(ENV) >> STAGED_PATH
        1 * cache.runCommand({ it.contains('mamba env create') && it.contains('--yes') && it.contains('--file') && it.contains(STAGED_PATH.toString()) }) >> null
        result == PREFIX
    }

    def 'should create a conda environment using micromamba and remote yaml file' () {
        given:
        def ENV = 'http://foo.com/some/env.yml'
        def PREFIX = Files.createTempDirectory('foo')
        def STAGED_PATH = Paths.get('/staged/env.yml')
        def cache = Spy(new CondaCache(useMicromamba: true))

        when:
        // the prefix directory exists ==> no mamba command is executed
        def result = cache.createLocalCondaEnv(ENV, PREFIX)
        then:
        0 * cache.isYamlFilePath(ENV)
        0 * cache.runCommand(_)
        result == PREFIX

        when:
        PREFIX.deleteDir()
        result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        // isYamlFilePath returns true for .yml files, so we enter YAML branch
        1 * cache.isYamlFilePath(ENV)
        // getLocalFilePath is called, which internally calls isRemoteFile and stageRemoteFile
        1 * cache.getLocalFilePath(ENV) >> STAGED_PATH
        1 * cache.runCommand({ it.contains('micromamba env create') && it.contains('--yes') && it.contains('--file') && it.contains(STAGED_PATH.toString()) }) >> null
        result == PREFIX
    }

    def 'should create a conda environment using mamba and remote lock file' () {
        given:
        // Use a URL without .yml extension so isYamlFilePath returns false
        def ENV = 'http://foo.com/some/condalock'
        def PREFIX = Files.createTempDirectory('foo')
        def STAGED_PATH = Paths.get('/staged/condalock')
        def cache = Spy(new CondaCache(useMamba: true))

        when:
        PREFIX.deleteDir()
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)  // returns false
        1 * cache.isRemoteFile(ENV)    // returns true, enters remote branch
        1 * cache.getLocalFilePath(ENV) >> STAGED_PATH
        1 * cache.isLockFilePath(STAGED_PATH) >> true
        1 * cache.runCommand({ it.contains('mamba create') && it.contains('--file') && it.contains(STAGED_PATH.toString()) }) >> null
        result == PREFIX
    }

    def 'should create a conda environment using micromamba and remote lock file' () {
        given:
        // Use a URL without .yml extension so isYamlFilePath returns false
        def ENV = 'http://foo.com/some/condalock'
        def PREFIX = Files.createTempDirectory('foo')
        def STAGED_PATH = Paths.get('/staged/condalock')
        def cache = Spy(new CondaCache(useMicromamba: true))

        when:
        PREFIX.deleteDir()
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)  // returns false
        1 * cache.isRemoteFile(ENV)    // returns true, enters remote branch
        1 * cache.getLocalFilePath(ENV) >> STAGED_PATH
        1 * cache.isLockFilePath(STAGED_PATH) >> true
        1 * cache.runCommand({ it.contains('micromamba create') && it.contains('--file') && it.contains(STAGED_PATH.toString()) }) >> null
        result == PREFIX
    }

    def 'should create conda env with options' () {
        given:
        def ENV = 'bwa=1.1.1'
        def PREFIX = Paths.get('/foo/bar')
        and:
        def cache = Spy(new CondaCache(createOptions: '--this --that'))

        when:
        def result = cache.createLocalCondaEnv0(ENV,PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.isTextFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand( "conda create --this --that --yes --quiet --prefix $PREFIX $ENV" ) >> null
        result == PREFIX
    }

    def 'should create conda env with options - using mamba' () {
        given:
        def ENV = 'bwa=1.1.1'
        def PREFIX = Paths.get('/foo/bar')
        and:
        def cache = Spy(new CondaCache(useMamba: true, createOptions: '--this --that'))

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.isTextFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand("mamba create --this --that --yes --quiet --prefix $PREFIX $ENV") >> null
        result == PREFIX
    }

    def 'should create conda env with options - using micromamba' () {
        given:
        def ENV = 'bwa=1.1.1'
        def PREFIX = Paths.get('/foo/bar')
        and:
        def cache = Spy(new CondaCache(useMicromamba: true, createOptions: '--this --that'))

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.isTextFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand("micromamba create --this --that --yes --quiet --prefix $PREFIX $ENV") >> null
        result == PREFIX
    }

    def 'should create conda env with channels' () {
        given:
        def ENV = 'bwa=1.1.1'
        def PREFIX = Paths.get('/foo/bar')
        and:
        def cache = Spy(new CondaCache(new CondaConfig([channels:['bioconda','defaults']], [:])))

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.isTextFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand("conda create --yes --quiet --prefix /foo/bar -c bioconda -c defaults bwa=1.1.1") >> null
        result == PREFIX
    }

    def 'should create a conda env with a yaml file' () {

        given:
        def ENV = 'foo.yml'
        def PREFIX = Paths.get('/conda/envs/my-env')
        def cache = Spy(CondaCache)

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        0 * cache.isTextFilePath(ENV)
        1 * cache.makeAbsolute(ENV) >> Paths.get('/usr/base').resolve(ENV)
        1 * cache.runCommand( "conda env create --prefix $PREFIX --file /usr/base/foo.yml" ) >> null
        result == PREFIX

    }

    def 'should create a conda env with a yaml file - using micromamba' () {

        given:
        def ENV = 'foo.yml'
        def PREFIX = Paths.get('/conda/envs/my-env')
        def cache = Spy(new CondaCache(useMicromamba: true))

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        0 * cache.isTextFilePath(ENV)
        1 * cache.makeAbsolute(ENV) >> Paths.get('/usr/base').resolve(ENV)
        1 * cache.runCommand( "micromamba env create --yes --prefix $PREFIX --file /usr/base/foo.yml" ) >> null
        result == PREFIX

    }

    def 'should create a conda env with a text file' () {

        given:
        def ENV = 'foo.txt'
        def PREFIX = Paths.get('/conda/envs/my-env')
        and:
        def cache = Spy(new CondaCache(createOptions: '--this --that'))

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.isTextFilePath(ENV)
        1 * cache.makeAbsolute(ENV) >> Paths.get('/usr/base').resolve(ENV)
        1 * cache.runCommand( "conda create --this --that --yes --quiet --prefix $PREFIX --file /usr/base/foo.txt" ) >> null
        result == PREFIX

    }

    def 'should create a conda env with a text file - using micromamba' () {

        given:
        def ENV = 'foo.txt'
        def PREFIX = Paths.get('/conda/envs/my-env')
        and:
        def cache = Spy(new CondaCache(useMicromamba: true, createOptions: '--this --that'))

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.isTextFilePath(ENV)
        1 * cache.makeAbsolute(ENV) >> Paths.get('/usr/base').resolve(ENV)
        1 * cache.runCommand( "micromamba create --this --that --yes --quiet --prefix $PREFIX --file /usr/base/foo.txt" ) >> null
        result == PREFIX

    }

    def 'should create a conda env with a lock file with .lock extension' () {

        given:
        def folder = Files.createTempDirectory('test')
        def ENV = folder.resolve('env.lock').toString()
        folder.resolve('env.lock').text = '''
                @EXPLICIT
                https://conda.anaconda.org/conda-forge/linux-64/package-1.0.tar.bz2
                '''
                .stripIndent(true)
        def PREFIX = Paths.get('/conda/envs/my-env')
        and:
        def cache = Spy(new CondaCache(createOptions: '--this --that'))

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.isRemoteFile(ENV)
        1 * cache.isTextFilePath(ENV)
        1 * cache.isLockFilePath(_) >> true
        1 * cache.runCommand({ it.contains('--file') && it.contains('env.lock') }) >> null
        result == PREFIX

        cleanup:
        folder?.deleteDir()
    }

    def 'should create a conda env with a lock file with no extension' () {

        given:
        def folder = Files.createTempDirectory('test')
        def ENV = folder.resolve('condalock').toString()
        folder.resolve('condalock').text = '''
                @EXPLICIT
                https://conda.anaconda.org/conda-forge/linux-64/package-1.0.tar.bz2
                '''
                .stripIndent(true)
        def PREFIX = Paths.get('/conda/envs/my-env')
        and:
        def cache = Spy(new CondaCache())

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.isRemoteFile(ENV)
        1 * cache.isTextFilePath(ENV)
        1 * cache.isLockFilePath(_) >> true
        1 * cache.runCommand({ it.contains('--file') && it.contains('condalock') }) >> null
        result == PREFIX

        cleanup:
        folder?.deleteDir()
    }

    def 'should get options from the config' () {

        when:
        def config = new CondaConfig([:], [:])
        def cache = new CondaCache(config)
        then:
        cache.createTimeout.minutes == 20
        cache.createOptions == null
        cache.configCacheDir0 == null
        !cache.@useMamba
        !cache.@useMicromamba
        cache.binaryName == "conda"

        when:
        config = new CondaConfig([createTimeout: '5 min', createOptions: '--foo --bar', cacheDir: '/conda/cache', useMamba: true], [:])
        cache = new CondaCache(config)
        then:
        cache.createTimeout.minutes == 5
        cache.createOptions == '--foo --bar'
        cache.configCacheDir0 == Paths.get('/conda/cache')
        cache.@useMamba
        !cache.@useMicromamba
        cache.binaryName == "mamba"

        when:
        config = new CondaConfig([createTimeout: '5 min', createOptions: '--foo --bar', cacheDir: '/conda/cache', useMicromamba: true], [:])
        cache = new CondaCache(config)
        then:
        cache.createTimeout.minutes == 5
        cache.createOptions == '--foo --bar'
        cache.configCacheDir0 == Paths.get('/conda/cache')
        !cache.@useMamba
        cache.@useMicromamba
        cache.binaryName == "micromamba"
    }

    def 'should define cache dir from config' () {

        given:
        def folder = Files.createTempDirectory('test'); folder.deleteDir()
        def config = new CondaConfig([cacheDir: folder.toString()], [:])
        CondaCache cache = Spy(CondaCache, constructorArgs: [config])

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
        def folder = Paths.get('.test-conda-cache-' + Math.random())
        def config = new CondaConfig([cacheDir: folder.toString()], [:])
        CondaCache cache = Spy(CondaCache, constructorArgs: [config])

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
        def config = new CondaConfig()
        CondaCache cache = Spy(CondaCache, constructorArgs: [config])

        when:
        def result = cache.getCacheDir()
        then:
        2 * cache.getEnv() >> [NXF_CONDA_CACHEDIR: folder.toString()]
        0 * cache.getSessionWorkDir()
        result == folder
        result.exists()

        cleanup:
        folder?.deleteDir()
    }

    def 'should define cache dir from session workdir' () {

        given:
        def folder = Files.createTempDirectory('test');
        def cache = Spy(CondaCache)

        when:
        def result = cache.getCacheDir()
        then:
        1 * cache.getSessionWorkDir() >> folder
        result == folder.resolve('conda')
        result.exists()

        cleanup:
        folder?.deleteDir()
    }
}
