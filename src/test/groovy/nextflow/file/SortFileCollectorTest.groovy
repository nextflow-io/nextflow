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
import java.nio.file.Files

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SortFileCollectorTest extends Specification {

    def testAppendString() {

        given:
        def folder = Files.createTempDirectory('test')

        when:
        def collector = new SortFileCollector()
        collector.append('alpha', 'BBB')
        collector.append('delta', '222')
        collector.append('alpha', 'AAA')
        collector.append('delta', '111')
        collector.moveFiles(folder)
        /*
         * No sorting has been specified
         * Entries are appended in the order they have been added
         */
        then:
        folder.list().size() == 2
        folder.resolve('alpha').text == 'BBBAAA'
        folder.resolve('delta').text == '222111'

        cleanup:
        collector?.close()
        folder?.deleteDir()

    }

    def testAppendStringWithSort() {

        given:
        def folder = Files.createTempDirectory('test')

        when:
        def collector = new SortFileCollector()
        collector.sort = { it -> it }
        collector.append('alpha', 'BBB')
        collector.append('delta', '222')
        collector.append('alpha', 'AAA')
        collector.append('delta', '111')
        collector.moveFiles(folder)
        then:
        folder.list().size() == 2
        folder.resolve('alpha').text == 'AAABBB'
        folder.resolve('delta').text == '111222'

        cleanup:
        collector?.close()
        folder?.deleteDir()

    }


    def testAppendStringWithSortSeedAndNewLine() {

        given:
        def folder = Files.createTempDirectory('test')

        when:
        def collector = new SortFileCollector()
        collector.sort = { it -> it }
        collector.newLine = true
        collector.seed = [alpha: '000', delta: '111', gamma: '222']
        collector.append('alpha', 'BBB')
        collector.append('delta', 'qqq')
        collector.append('alpha', 'ZZZ')
        collector.append('delta', 'ttt')
        collector.append('gamma', 'yyy')
        collector.append('gamma', 'zzz')
        collector.append('gamma', 'xxx')
        collector.append('delta', 'ppp')
        collector.append('alpha', 'AAA')
        collector.moveFiles(folder)
        then:
        folder.list().size() == 3
        folder.resolve('alpha').text == '000\nAAA\nBBB\nZZZ\n'
        folder.resolve('delta').text == '111\nppp\nqqq\nttt\n'
        folder.resolve('gamma').text == '222\nxxx\nyyy\nzzz\n'

        cleanup:
        collector?.close()
        folder?.deleteDir()

    }
}
