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
class ShifterBuilderTest extends Specification {

    def 'test shifter env'() {

        given:
        def builder = new ShifterBuilder('x')

        expect:
        builder.makeEnv('X=1').toString() == 'X=1'
        builder.makeEnv([VAR_X:1, VAR_Y: 2]).toString() == 'VAR_X="1" VAR_Y="2"'
    }

    def 'should build the shifter run command' () {

        expect:
        new ShifterBuilder('busybox')
                .build()
                .@runCommand == 'shifter --image busybox'

        new ShifterBuilder('busybox')
                .params(verbose: true)
                .build()
                .@runCommand == 'shifter --verbose --image busybox'

        new ShifterBuilder('fedora')
                .addEnv([VAR_X:1, VAR_Y:2])
                .addEnv("VAR_Z=3")
                .build()
                .@runCommand == 'VAR_X="1" VAR_Y="2" VAR_Z=3 shifter --image fedora'

    }

    def 'should get run command line' () {

        when:
        def cli = new ShifterBuilder('ubuntu:14').build().getRunCommand()
        then:
        cli ==  '''
                shifter_pull ubuntu:14
                shifter --image ubuntu:14
                '''
                .stripIndent().trim()

        when:
        cli = new ShifterBuilder('ubuntu:14').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  '''
                shifter_pull ubuntu:14
                shifter --image ubuntu:14 bwa --this --that file.fasta
                '''
                .stripIndent().trim()

        when:
        cli = new ShifterBuilder('ubuntu:14').params(entry:'/bin/bash').build().getRunCommand('bwa --this --that file.fasta')
        then:
        cli ==  '''
                shifter_pull ubuntu:14
                shifter --image ubuntu:14 /bin/bash -c "bwa --this --that file.fasta"
                '''
                .stripIndent().trim()

    }



}
