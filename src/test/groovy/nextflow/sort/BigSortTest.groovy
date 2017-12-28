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

package nextflow.sort

import java.nio.file.Paths

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BigSortTest extends Specification {

    def testProperties () {

        given:
        def testComp = new Comparator() {
            @Override
            int compare(Object o1, Object o2) {
                return 0
            }
        }

        when:
        def sort = [:] as BigSort
        sort.comparator(testComp)
        sort.sliceMaxItems(1_000)
        sort.sliceMaxSize(2_000)
        sort.tempDir(Paths.get('/some/path'))
        sort.deleteTempFilesOnClose(false)

        then:
        sort.comparator == testComp
        sort.sliceMaxItems == 1_000
        sort.sliceMaxSize == 2_000
        sort.tempDir == Paths.get('/some/path')
        sort.deleteTempFilesOnClose == false
    }


    def testSorting() {

        given:
        def sort = new BigSort<String>() {

            def Map<Long,String> map = new HashMap<>()

            @Override
            protected int put(long key, String value) {
                map.put(key,value)
                return 8 + value.length()
            }

            @Override
            protected String get(long key) {
                return map.get(key)
            }
        }

        sort.sliceMaxItems(3)
        sort.create()

        when:
        sort.add('444')
        sort.add('333')
        sort.add('999')
        sort.add('666')
        sort.add('777')
        sort.add('222')
        sort.add('111')
        sort.add('555')
        sort.add('000')
        sort.add('888')

        def result = []
        sort.sort { result.add(it) }

        then:
        // first slice is at 0 (implicit), second starts at 3, third at 6, and last one with at position 9
        sort.slices == [3, 6, 9]
        result.size() == 10
        result == ['000','111','222','333','444','555','666','777','888','999']

        cleanup:
        sort?.close()

    }

}
