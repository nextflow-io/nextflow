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

package nextflow.util

import spock.lang.Specification
import spock.lang.Timeout

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class NameGeneratorTest extends Specification {

    def 'should return a random name' () {
        when:
        def (adj, name)= NameGenerator.next().tokenize('_')
        then:
        NameGenerator.ADJECTIVES.contains(adj)
        NameGenerator.NAMES.contains(name)

    }

    @Timeout(1)
    def 'should not generate a random name except the specified one' () {
        when:
        def name = NameGenerator.next('evil_pike')
        then:
        name
        name != 'evil_pike'
    }
}
