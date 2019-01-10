/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
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
