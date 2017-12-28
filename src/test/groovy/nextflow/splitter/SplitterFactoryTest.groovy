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

package nextflow.splitter

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SplitterFactoryTest extends Specification {

    def testCreateSplitter() {

        expect:
        SplitterFactory.create('text') instanceof TextSplitter
        SplitterFactory.create('fasta') instanceof FastaSplitter

    }

    def testArgsToOptions() {

        given:
        def closure = { -> 1 }

        expect:
        SplitterFactory.argsToOpt( [ ] as Object[] ) == [:]
        SplitterFactory.argsToOpt( [ closure ] as Object[] ) == [ each: closure ]
        SplitterFactory.argsToOpt( [ [x:1, y:2] ] as Object[] ) == [x:1, y:2]
        SplitterFactory.argsToOpt( [ [x:1, y:2], closure ] as Object[] ) == [x:1, y:2, each: closure]
    }

}
