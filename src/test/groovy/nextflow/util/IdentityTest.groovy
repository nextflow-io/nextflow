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

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class IdentityTest extends Specification{

    def 'should validate hash code' () {

        expect:
        new HexIdentity(1).hashCode() == new HexIdentity(1).hashCode()
        new HexIdentity(100).hashCode() == new HexIdentity(100).hashCode()
        new HexIdentity(100).hashCode() != new HexIdentity(101).hashCode()
    }


    def "should validate equals methods" () {

        expect:
        new HexIdentity(1) == new HexIdentity(1)
        new HexIdentity(100) == new HexIdentity(100)
        new HexIdentity(100) != new HexIdentity(101)

    }

    def "should validate toString method" () {

        expect:
        new HexIdentity(1).toString() == "0x1"
        new HexIdentity(10).toString() == "0xa"
        new HexIdentity(100).toString() == "0x64"

    }

}
