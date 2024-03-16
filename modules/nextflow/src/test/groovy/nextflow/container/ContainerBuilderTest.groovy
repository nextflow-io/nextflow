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

import java.nio.file.Path

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerBuilderTest extends Specification {

    def 'should return mount flags'() {

        given:
        def builder = Spy(ContainerBuilder)

        expect:
        builder.mountFlags(false) == ''
        builder.mountFlags(true) == ':ro'

    }

    def 'should make env var' () {
        given:
        StringBuilder result
        def builder = Spy(ContainerBuilder)

        when:
        result = builder.makeEnv([FOO: 'x', BAR: 'y'])
        then:
        result.toString() == '-e "FOO=x" -e "BAR=y"'

        when:
        result = builder.makeEnv('FOO=hello')
        then:
        result.toString() == '-e "FOO=hello"'
        
        when:
        result = builder.makeEnv( 'FOO' )
        then:
        result.toString() == '-e "FOO"'

        when:
        builder.makeEnv( 1 )
        then:
        thrown(IllegalArgumentException)

    }

    @Unroll
    def 'should create builder for given engine' () {
        given:
        def IMAGE = 'foo:latest'

        when:
        def builder = ContainerBuilder.create(ENGINE,IMAGE)
        then:
        builder.class == CLAZZ
        builder.getImage() == IMAGE

        where:
        ENGINE              | CLAZZ
        'docker'            | DockerBuilder
        'podman'            | PodmanBuilder
        'singularity'       | SingularityBuilder
        'apptainer'         | ApptainerBuilder
        'sarus'             | SarusBuilder
        'shifter'           | ShifterBuilder
        'charliecloud'      | CharliecloudBuilder
        'udocker'           | UdockerBuilder

    }

    def 'should throw illegal arg' () {

        when:
        ContainerBuilder.create('foo','image:any')

        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Unknown container engine: foo'

    }

    def 'should add mount'() {
        given:
        def builder = Spy(ContainerBuilder)

        when:
        builder.addMount(null)
        builder.addMount(Path.of('/this'))
        builder.addMount(Path.of('/that'))
        then:
        builder.@mounts == [Path.of('/this'), Path.of('/that')]
    }

    def 'should add mounts'() {
        given:
        def builder = Spy(ContainerBuilder)

        when:
        builder.addMounts(null)
        builder.addMounts([Path.of('/this'), Path.of('/that')])
        builder.addMounts([Path.of('/another'), Path.of('/one')])
        then:
        builder.@mounts == [Path.of('/this'), Path.of('/that'), Path.of('/another'), Path.of('/one')]

    }
}
