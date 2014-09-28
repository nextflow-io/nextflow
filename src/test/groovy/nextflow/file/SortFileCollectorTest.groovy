/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

package nextflow.file

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SortFileCollectorTest extends Specification {

    def testPut() {

        given:
        def collector = new SortFileCollector<Integer,String>()

        when:
        collector.put( 10, 'Hello' )
        then:
        collector.getStore().get(0) == 'Hello'
        collector.getIndex().get(0).key == 10
        collector.getIndex().get(0).pos == 0
        collector.getIndex().get(0).left == -1
        collector.getIndex().get(0).right == -1
        collector.getCount() == 1

        when:
        collector.put( 2, 'Hola' )
        then:
        collector.getStore().get(0) == 'Hello'
        collector.getIndex().get(0).key == 10
        collector.getIndex().get(0).pos == 0
        collector.getIndex().get(0).left == 1
        collector.getIndex().get(0).right == -1

        collector.getStore().get(1) == 'Hola'
        collector.getIndex().get(1).key == 2
        collector.getIndex().get(1).pos == 1
        collector.getIndex().get(1).left == -1
        collector.getIndex().get(1).right == -1

        collector.getCount() == 2

        when:
        collector.put( 5, 'Ciao' )
        then:
        collector.getStore().get(0) == 'Hello'
        collector.getIndex().get(0).key == 10
        collector.getIndex().get(0).pos == 0
        collector.getIndex().get(0).left == 1
        collector.getIndex().get(0).right == -1

        collector.getStore().get(1) == 'Hola'
        collector.getIndex().get(1).key == 2
        collector.getIndex().get(1).pos == 1
        collector.getIndex().get(1).left == -1
        collector.getIndex().get(1).right == 2

        collector.getStore().get(2) == 'Ciao'
        collector.getIndex().get(2).key == 5
        collector.getIndex().get(2).pos == 2
        collector.getIndex().get(2).left == -1
        collector.getIndex().get(2).right == -1

        collector.getCount() == 3

        // overwrite value for '2'
        when:
        collector.put( 2, 'Bye' )
        then:
        collector.getStore().get(0) == 'Hello'
        collector.getIndex().get(0).key == 10
        collector.getIndex().get(0).pos == 0
        collector.getIndex().get(0).left == 1
        collector.getIndex().get(0).right == -1

        collector.getStore().get(1) == 'Bye'
        collector.getIndex().get(1).key == 2
        collector.getIndex().get(1).pos == 1
        collector.getIndex().get(1).left == -1
        collector.getIndex().get(1).right == 2

        collector.getStore().get(2) == 'Ciao'
        collector.getIndex().get(2).key == 5
        collector.getIndex().get(2).pos == 2
        collector.getIndex().get(2).left == -1
        collector.getIndex().get(2).right == -1

        collector.getCount() == 3

    }

    def testGet() {

        given:
        def collector = new SortFileCollector<Integer,String>()

        when:
        collector.put( 9, 'Hello' )
        collector.put( 3, 'Ciao' )
        collector.put( 5, 'Hola' )
        collector.put( 1, 'Bye' )
        collector.put( 7, 'Bonjour' )

        then:
        collector.get( 1 ) == 'Bye'
        collector.get( 3 ) == 'Ciao'
        collector.get( 5 ) == 'Hola'
        collector.get( 7 ) == 'Bonjour'
        collector.get( 9 ) == 'Hello'
        collector.get( 2 ) == null

    }

    def testContainsKey() {

        given:
        def collector = new SortFileCollector<Integer,String>()

        when:
        collector.put( 9, 'Hello' )
        collector.put( 3, 'Ciao' )
        collector.put( 5, 'Hola' )
        collector.put( 1, 'Bye' )
        collector.put( 7, 'Bonjour' )

        then:
        collector.containsKey(1)
        collector.containsKey(3)
        collector.containsKey(5)
        collector.containsKey(7)
        !collector.containsKey(99)

    }


    def testIterator() {

        given:
        def collector = new SortFileCollector<Integer,String>()
        collector.put( 9, 'Hello' )
        collector.put( 3, 'Ciao' )
        collector.put( 5, 'Hola' )
        collector.put( 1, 'Bye' )
        collector.put( 7, 'Bonjour' )

        when:
        def itr = collector.iterator()

        then:
        itr.hasNext()
        itr.next() ==  1
        itr.hasNext()
        itr.next() ==  3
        itr.hasNext()
        itr.next() ==  5
        itr.hasNext()
        itr.next() ==  7
        itr.hasNext()
        itr.next() ==  9
        !itr.hasNext()
        !itr.hasNext()
        itr.next() == null
        itr.next() == null


    }

}
