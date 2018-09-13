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

package nextflow.container

import spock.lang.Specification

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
        result.toString() == '${FOO:+-e "FOO=$FOO"}'

        when:
        builder.makeEnv( 1 )
        then:
        thrown(IllegalArgumentException)

    }

}
