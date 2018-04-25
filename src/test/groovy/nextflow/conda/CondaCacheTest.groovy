/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
 */

package nextflow.conda

import java.nio.file.Files
import java.nio.file.Paths

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CondaCacheTest extends Specification {

    def 'should env file' () {

        given:
        def cache = new CondaCache()
        
        expect:
        !cache.isEnvFile('foo=1.0')
        cache.isEnvFile('env.yml')
        cache.isEnvFile('env.yaml')
    }


    def 'should create conda env prefix path for a string env' () {

        given:
        def ENV = 'bwa=1.7.2'
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')

        when:
        def prefix = cache.condaPrefixPath(ENV)
        then:
        1 * cache.isEnvFile(ENV)
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == '/conda/envs/env-eaeb133f4ca62c95e9c0eec7ef8d553b'
    }

    def 'should create conda env prefix path for a env file' () {

        given:
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')
        def ENV = Files.createTempFile('test','.yml')
        ENV.text = '''
            channels:
              - bioconda
              - defaults
            dependencies:
              # Default bismark
              - star=2.5.4a
              - bwa=0.7.15        
            '''
            .stripIndent()

        when:
        def prefix = cache.condaPrefixPath(ENV.toString())
        then:
        1 * cache.isEnvFile(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == '/conda/envs/env-9416240708c49c4e627414b46a743664'

    }

    def 'should create conda env prefix path for a env file with name' () {

        given:
        def cache = Spy(CondaCache)
        def BASE = Paths.get('/conda/envs')
        def ENV = Files.createTempFile('test','.yml')
        ENV.text = '''  
            name: my-env-1.1
            channels:
              - bioconda
              - defaults
            dependencies:
              # Default bismark
              - star=2.5.4a
              - bwa=0.7.15        
            '''
                .stripIndent()

        when:
        def prefix = cache.condaPrefixPath(ENV.toString())
        then:
        1 * cache.isEnvFile(ENV.toString())
        1 * cache.getCacheDir() >> BASE
        prefix.toString() == '/conda/envs/my-env-1.1-e7fafe40ca966397a2c0d9bed7181aa7'

    }


    def 'should create a conda environment' () {

        given:
        def ENV = 'bwa=1.1.1'
        def PREFIX = Files.createTempDirectory('foo')
        def cache = Spy(CondaCache)

        when:
        // the prefix directory exists ==> no conda command is executed
        def result = cache.createLocalCondaEnv0(ENV,PREFIX)
        then:
        0 * cache.isEnvFile(ENV)
        0 * cache.runCommand(_)
        result == PREFIX

        when:
        PREFIX.deleteDir()
        result = cache.createLocalCondaEnv0(ENV,PREFIX)
        then:
        1 * cache.isEnvFile(ENV)
        0 * cache.makeAbsolute(_)
        1 * cache.runCommand( "conda create --mkdir --yes --quiet --prefix $PREFIX $ENV" ) >> null
        result == PREFIX
        
    }

    def 'should create a conda env with a env file' () {

        given:
        def ENV = 'foo.yml'
        def PREFIX = Paths.get('/conda/envs/my-env')
        def cache = Spy(CondaCache)

        when:
        def result = cache.createLocalCondaEnv0(ENV, PREFIX)
        then:
        1 * cache.isEnvFile(ENV)
        1 * cache.makeAbsolute(ENV) >> Paths.get('/usr/base').resolve(ENV)
        1 * cache.runCommand( "conda env create --prefix $PREFIX --file /usr/base/foo.yml" ) >> null
        result == PREFIX

    }
}
