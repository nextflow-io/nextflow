/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
class HistoryFileTest extends Specification {

    def 'test add and get and find' () {

        setup:
        new File('.nextflow.history').delete()

        when:
        true

        then:
        HistoryFile.history.retrieveLastUniqueId() == null

        when:
        def id1 = UUID.randomUUID()
        def id2 = UUID.randomUUID()
        def id3 = UUID.randomUUID()
        HistoryFile.history.write( id1, [1,2,3] )
        HistoryFile.history.write( id2, [1,2,3] )
        HistoryFile.history.write( id3, [1,2,3] )

        then:
        HistoryFile.history.retrieveLastUniqueId() == id3.toString()
        HistoryFile.history.findUniqueId( id1.toString() )
        HistoryFile.history.findUniqueId( id2.toString() )
        HistoryFile.history.findUniqueId( id3.toString() )
        !HistoryFile.history.findUniqueId( UUID.randomUUID().toString() )

        cleanup:
        HistoryFile.history.delete()

    }


}
