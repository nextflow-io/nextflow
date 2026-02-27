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

    def 'should detect explicit file by content' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cache = new CondaCache()

        // File with @EXPLICIT marker
        def explicitFile = folder.resolve('packages.txt')
        explicitFile.text = '''\
            # This file was created by conda
            @EXPLICIT
            https://repo.anaconda.com/pkgs/main/linux-64/bwa-0.7.17.tar.bz2
            '''.stripIndent()

        // File without @EXPLICIT marker
        def regularFile = folder.resolve('regular.txt')
        regularFile.text = '''\
            bwa=0.7.17
            samtools=1.9
            '''.stripIndent()

        // YAML file (should not match)
        def yamlFile = folder.resolve('env.yaml')
        yamlFile.text = '''\
            dependencies:
              - bwa=0.7.17
            '''.stripIndent()

        expect:
        // Non-existent file
        !cache.isExplicitFile('foo=1.0')
        !cache.isExplicitFile('nonexistent.txt')
        // String with newline
        !cache.isExplicitFile('foo.txt\nbar.txt')
        // File with @EXPLICIT marker
        cache.isExplicitFile(explicitFile.toString())
        // File without @EXPLICIT marker
        !cache.isExplicitFile(regularFile.toString())
        // YAML file
        !cache.isExplicitFile(yamlFile.toString())

        cleanup:
        folder?.deleteDir()
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

        when:
        def prefix = cache.condaPrefixPath(ENV)
        then:
        0 * cache.isYamlFilePath(ENV)
        1 * cache.isYamlUriPath(ENV)
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == '/conda/envs/env-12c863103deed9425ce8012323f948fc'
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

    def 'should create conda env prefix path for a conda explicit file' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')
        def ENV = folder.resolve('bar.txt')
        ENV.text = '''\
                # This file was created by conda
                @EXPLICIT
                https://repo.anaconda.com/pkgs/main/linux-64/star-2.5.4a.tar.bz2
                https://repo.anaconda.com/pkgs/main/linux-64/bwa-0.7.15.tar.bz2
                https://repo.anaconda.com/pkgs/main/linux-64/multiqc-1.2.3.tar.bz2
                '''.stripIndent()

        when:
        def prefix = cache.condaPrefixPath(ENV.toString())
        then:
        1 * cache.isYamlFilePath(ENV.toString())
        1 * cache.isExplicitFile(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == "/conda/envs/env-24d602a5eecba868858ab48a41e2c9bd"

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

    def 'should create a conda environment using mamba and remote lock file' () {
        given:
        def ENV = 'http://foo.com/some/file-lock.yml'
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
        1 * cache.runCommand("mamba env create --yes --prefix $PREFIX --file $ENV") >> null
        result == PREFIX
    }

    def 'should create a conda environment using micromamba and remote lock file' () {
        given:
        def ENV = 'http://foo.com/some/file-lock.yml'
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
        1 * cache.runCommand("micromamba env create --yes --prefix $PREFIX --file $ENV") >> null
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
        1 * cache.isExplicitFile(ENV)
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
        1 * cache.isExplicitFile(ENV)
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
        1 * cache.isExplicitFile(ENV)
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
        1 * cache.isExplicitFile(ENV)
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
        0 * cache.isExplicitFile(ENV)
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
        0 * cache.isExplicitFile(ENV)
        1 * cache.makeAbsolute(ENV) >> Paths.get('/usr/base').resolve(ENV)
        1 * cache.runCommand( "micromamba env create --yes --prefix $PREFIX --file /usr/base/foo.yml" ) >> null
        result == PREFIX

    }

    def 'should create a conda env with an explicit file' () {

        given:
        def folder = Files.createTempDirectory('test')
        def envFile = folder.resolve('explicit.txt')
        envFile.text = '''\
            # This file may be used to create an environment using:
            # $ conda create --name <env> --file <this file>
            # platform: linux-64
            @EXPLICIT
            https://conda.anaconda.org/conda-forge/linux-64/_libgcc_mutex-0.1-conda_forge.tar.bz2#d7c89558ba9fa0495403155b64376d81
            https://conda.anaconda.org/conda-forge/linux-64/libgomp-15.2.0-h767d61c_7.conda#f7b4d76975aac7e5d9e6ad13845f92fe
            https://conda.anaconda.org/conda-forge/linux-64/_openmp_mutex-4.5-2_gnu.tar.bz2#73aaf86a425cc6e73fcf236a5a46396d
            https://conda.anaconda.org/conda-forge/linux-64/libgcc-15.2.0-h767d61c_7.conda#c0374badb3a5d4b1372db28d19462c53
            '''.stripIndent()
        def ENV = envFile.toString()
        def PREFIX = Paths.get('/conda/envs/my-env')
        and:
        def cache = Spy(new CondaCache(createOptions: '--this --that'))

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.isExplicitFile(ENV)
        1 * cache.makeAbsolute(ENV) >> envFile.toAbsolutePath()
        1 * cache.runCommand( "conda create --this --that --yes --quiet --prefix $PREFIX --file ${envFile.toAbsolutePath()}" ) >> null
        result == PREFIX

        cleanup:
        folder?.deleteDir()
    }

    def 'should create a conda env with an explicit file - using micromamba' () {

        given:
        def folder = Files.createTempDirectory('test')
        def envFile = folder.resolve('explicit.txt')
        envFile.text = '''\
            # This file may be used to create an environment using:
            # $ conda create --name <env> --file <this file>
            # platform: linux-64
            @EXPLICIT
            https://conda.anaconda.org/conda-forge/linux-64/_libgcc_mutex-0.1-conda_forge.tar.bz2#d7c89558ba9fa0495403155b64376d81
            https://conda.anaconda.org/conda-forge/linux-64/libgomp-15.2.0-h767d61c_7.conda#f7b4d76975aac7e5d9e6ad13845f92fe
            https://conda.anaconda.org/conda-forge/linux-64/_openmp_mutex-4.5-2_gnu.tar.bz2#73aaf86a425cc6e73fcf236a5a46396d
            '''.stripIndent()
        def ENV = envFile.toString()
        def PREFIX = Paths.get('/conda/envs/my-env')
        and:
        def cache = Spy(new CondaCache(useMicromamba: true, createOptions: '--this --that'))

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.isExplicitFile(ENV)
        1 * cache.makeAbsolute(ENV) >> envFile.toAbsolutePath()
        1 * cache.runCommand( "micromamba create --this --that --yes --quiet --prefix $PREFIX --file ${envFile.toAbsolutePath()}" ) >> null
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
