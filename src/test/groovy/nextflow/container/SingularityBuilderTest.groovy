/*
 * Copyright (c) 2013-2016, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2016, Paolo Di Tommaso and the respective authors.
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

import java.nio.file.Paths

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SingularityBuilderTest extends Specification {

    def 'should get the exec command line' () {

        expect:
        new SingularityBuilder('busybox')
                .build()
                .runCommand == 'env - PATH="$PATH" singularity exec busybox'

        new SingularityBuilder('busybox')
                .params(entry: '/bin/bash')
                .build()
                .runCommand == 'env - PATH="$PATH" singularity exec busybox /bin/bash'

        new SingularityBuilder('busybox')
                .params(engineOptions: '-q -v')
                .build()
                .runCommand == 'env - PATH="$PATH" singularity -q -v exec busybox'

        new SingularityBuilder('busybox')
                .params(runOptions: '--contain --writable')
                .build()
                .runCommand == 'env - PATH="$PATH" singularity exec --contain --writable busybox'

    }

    def 'should return export variables' () {

        expect:
        new SingularityBuilder('busybox')
                .getEnvExports() == ''

        new SingularityBuilder('busybox')
                .addEnv('X=1')
                .addEnv(ALPHA:'aaa', BETA: 'bbb')
                .addEnv( new File('/bash/env.txt') )
                .getEnvExports() == 'export X=1; export ALPHA="aaa";  export BETA="bbb"; export BASH_ENV="/bash/env.txt"; '
    }

    def 'should normalise container path' () {
        expect:
        SingularityBuilder.normalizeImageName('/abs/path/bar.img') == '/abs/path/bar.img'
        SingularityBuilder.normalizeImageName('foo.img') == Paths.get('foo.img').toAbsolutePath().toString()
    }
}
