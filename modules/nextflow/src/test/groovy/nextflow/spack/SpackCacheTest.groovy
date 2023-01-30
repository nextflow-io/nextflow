/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.spack

import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Specification
/**
 *
 * @author Marco De La Pierre <marco.delapierre@gmail.com>
 */
class SpackCacheTest extends Specification {

    def 'should env file' () {

        given:
        def cache = new SpackCache()

        expect:
        !cache.isYamlFilePath('foo@1.0')
        cache.isYamlFilePath('env.yaml')
    }


    def 'should create spack env prefix path for a string env' () {

        given:
        def ENV = 'bwa@1.7.2'
        def cache = Spy(SpackCache)
        def BASE = Paths.get('/spack/envs')

        when:
        def prefix = cache.spackPrefixPath(ENV)
        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == '/spack/envs/env-dce1faa0a04e92386f05195016ca8440'
    }


    def 'should create spack env prefix path for a yaml env file' () {

        given:
        def folder = Files.createTempDirectory('test')
        def cache = Spy(SpackCache)
        def BASE = Paths.get('/spack/envs')
        def ENV = folder.resolve('foo.yaml')
        ENV.text = '''
            spack:
              specs: [star@2.5.4a, bwa@0.7.15]

              view: true
              concretizer:
                unify: true
            '''
            .stripIndent(true)  // https://issues.apache.org/jira/browse/GROOVY-9423

        when:
        def prefix = cache.spackPrefixPath(ENV.toString())
        then:
        1 * cache.isYamlFilePath(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == '/spack/envs/foo-ce5aa4c05c80124f5f5b75d4df55eeb5'

        cleanup:
        folder?.deleteDir()

    }

    def 'should return a spack prefix directory' () {

        given:
        def cache = Spy(SpackCache)
        def folder = Files.createTempDirectory('test')
        def ENV = folder.toString()

        when:
        def prefix = cache.spackPrefixPath(ENV)
        then:
        1 * cache.isYamlFilePath(ENV)
        0 * cache.getCacheDir()
        prefix.toString() == folder.toString()

        cleanup:
        folder?.deleteDir()

    }


    def 'should create a spack environment' () {

        given:
        def ENV = 'bwa@1.1.1'
        def PREFIX = Files.createTempDirectory('foo')
        def cache = Spy(SpackCache)

        when:
        // the prefix directory exists ==> no spack command is executed
        def result = cache.createLocalSpackEnv(ENV)

        then:
        1 * cache.spackPrefixPath(ENV) >> PREFIX
        0 * cache.isYamlFilePath(ENV)
        1 * cache.runCommand( "spack env activate $PREFIX ; spack install -y ; spack env deactivate" ) >> null
        result == PREFIX

        when:
        PREFIX.deleteDir()
        result = cache.createLocalSpackEnv0(ENV,PREFIX)
        then:
        1 * cache.isYamlFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand( "spack env create -d $PREFIX ; spack env activate $PREFIX ; spack add $ENV ; spack concretize -f ; spack install -y ; spack env deactivate" ) >> null
        result == PREFIX

    }

    def 'should create spack env with options' () {
        given:
        def ENV = 'bwa@1.1.1'
        def PREFIX = Paths.get('/foo/bar')
        and:
        def cache = Spy(new SpackCache([parallelBuilds: 2, noChecksum: true]))

        when:
        def result = cache.createLocalSpackEnv0(ENV,PREFIX)

        then:
        1 * cache.isYamlFilePath(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand( "spack env create -d $PREFIX ; spack env activate $PREFIX ; spack add $ENV ; spack concretize -f ; spack install -n -j 2 -y ; spack env deactivate" ) >> null
        result == PREFIX
    }

    def 'should create a spack env with a yaml file' () {

        given:
        def ENV = 'foo.yaml'
        def PREFIX = Paths.get('/spack/envs/my-env')
        def cache = Spy(SpackCache)

        when:
        def result = cache.createLocalSpackEnv0(ENV, PREFIX)

        then:
        1 * cache.isYamlFilePath(ENV)
        1 * cache.makeAbsolute(ENV) >> Paths.get('/usr/base').resolve(ENV)
        1 * cache.runCommand( "spack env create -d $PREFIX /usr/base/$ENV ; spack env activate $PREFIX ; spack concretize -f ; spack install -y ; spack env deactivate" ) >> null
        result == PREFIX

    }

    def 'should get options from the config' () {

        when:
        def cache = new SpackCache(new SpackConfig())
        then:
        cache.createTimeout.minutes == 60
        cache.configCacheDir0 == null
        !cache.@noChecksum
        cache.parallelBuilds == null

        when:
        cache = new SpackCache(new SpackConfig(createTimeout: '5 min', cacheDir: '/spack/cache', noChecksum: true, parallelBuilds: 2))
        then:
        cache.createTimeout.minutes == 5
        cache.configCacheDir0 == Paths.get('/spack/cache')
        cache.@noChecksum
        cache.parallelBuilds == 2
    }

    def 'should define cache dir from config' () {

        given:
        def folder = Files.createTempDirectory('test'); folder.deleteDir()
        def config = new SpackConfig(cacheDir: folder.toString())
        SpackCache cache = Spy(SpackCache, constructorArgs: [config])

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
        def folder = Paths.get('.test-spack-cache-' + Math.random())
        def config = new SpackConfig(cacheDir: folder.toString())
        SpackCache cache = Spy(SpackCache, constructorArgs: [config])

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
        def config = new SpackConfig()
        SpackCache cache = Spy(SpackCache, constructorArgs: [config])

        when:
        def result = cache.getCacheDir()
        then:
        2 * cache.getEnv() >> [NXF_SPACK_CACHEDIR: folder.toString()]
        0 * cache.getSessionWorkDir()
        result == folder
        result.exists()

        cleanup:
        folder?.deleteDir()
    }

    def 'should define cache dir from session workdir' () {

        given:
        def folder = Files.createTempDirectory('test');
        def cache = Spy(SpackCache)

        when:
        def result = cache.getCacheDir()
        then:
        1 * cache.getSessionWorkDir() >> folder
        result == folder.resolve('spack')
        result.exists()

        cleanup:
        folder?.deleteDir()
    }
}
