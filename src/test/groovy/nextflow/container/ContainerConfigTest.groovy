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
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ContainerConfigTest extends Specification {

    @Unroll
    def 'should return env whitelist for=#VAL' () {
        when:
        def cfg = new ContainerConfig(envWhitelist: VAL)
        then:
        cfg.getEnvWhitelist() == EXPECTED

        where:
        VAL         | EXPECTED
        null        | []
        ''          | []
        'FOO'       | ['FOO']
        'FOO,BAR'   | ['FOO','BAR']
        'A ,, B,C ' | ['A','B','C']
        ['X','Y']   | ['X','Y']

    }

}
