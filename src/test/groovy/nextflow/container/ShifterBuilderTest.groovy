/*
 * Copyright (c) 2013-2017, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2017, Paolo Di Tommaso and the respective authors.
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

import java.nio.file.Paths

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
        builder.makeEnv([VAR_X:1, VAR_Y: 2]).toString() == 'VAR_X=1 VAR_Y=2'
        builder.makeEnv( Paths.get('/some/file.env') ).toString() == 'BASH_ENV="/some/file.env"'
        builder.makeEnv( new File('/some/file.env') ).toString() == 'BASH_ENV="/some/file.env"'
    }

    def 'should build the shifter run command' () {

        expect:
        new ShifterBuilder('busybox')
                .build()
                .runCommand == 'shifter --image busybox'

        new ShifterBuilder('busybox')
                .params(verbose: true)
                .build()
                .runCommand == 'shifter --verbose --image busybox'

        new ShifterBuilder('ubuntu:latest')
                .params(entry: '/bin/bash')
                .build()
                .runCommand == 'shifter --image ubuntu:latest /bin/bash'

        new ShifterBuilder('ubuntu')
                .params(entry: '/bin/bash')
                .addEnv(Paths.get("/data/env_file"))
                .build()
                .runCommand == 'BASH_ENV="/data/env_file" shifter --image ubuntu /bin/bash'

        new ShifterBuilder('fedora')
                .params(entry: '/bin/bash')
                .addEnv([VAR_X:1, VAR_Y:2])
                .addEnv("VAR_Z=3")
                .build()
                .runCommand == 'VAR_X=1 VAR_Y=2 VAR_Z=3 shifter --image fedora /bin/bash'

    }


}
