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
class PathTrieTest extends Specification {

    def testAdd() {
        when:
        def tree = new PathTrie()
        tree.add('/home/data/work')
        tree.add('/db/blast/unitprot')
        then:
        tree.paths.size() == 2

        when:
        tree = new PathTrie()
        tree.add('/home/data/work')
        tree.add('/home/db/data')
        then:
        tree.paths.size() == 1
    }

    def testLongest() {

        when:
        def tree = new PathTrie()
        tree.add('/home/data/work')
        then:
        tree.longest().size() == 1
        tree.longest().get(0) == '/home/data/work'

        when:
        tree = new PathTrie()
        tree.add('/home/data/work')
        tree.add('/home/data/db')
        then:
        tree.longest().size() == 1
        tree.longest().get(0) == '/home/data'

        when:
        tree = new PathTrie()
        tree.add('/home/data/work')
        tree.add('/home/data/work/xx/file_x')
        tree.add('/home/data/work/yy/file_y')
        tree.add('/home/data/work/yy/zz/file_q')
        tree.add('/db/data/tutorial')
        tree.add('/usr/local/bin')
        then:
        tree.longest().size() == 3
        tree.longest().get(0) == '/home/data/work'
        tree.longest().get(1) == '/db/data/tutorial'
        tree.longest().get(2) == '/usr/local/bin'

        when:
        tree = new PathTrie()
        tree.add('a/b/c')
        tree.add('a/b/x')
        tree.add('p/q/z')
        then:
        tree.longest().size() == 2
        tree.longest().get(0) == 'a/b'
        tree.longest().get(1) == 'p/q/z'

    }

    def testLongest2 () {

        when:
        def tree = new PathTrie()
        tree.add('/home/data/')
        tree.add('/home/data/some/other/path')
        then:
        tree.longest().size() ==1
        tree.longest().get(0) == '/home/data'

        when:
        tree = new PathTrie()
        tree.add('/home/data/some/other/path')
        tree.add('/home/data')
        then:
        tree.longest().size() ==1
        tree.longest().get(0) == '/home/data'

        when:
        tree = new PathTrie()
        tree.add('/home/data/')
        tree.add('/home/data/some')
        tree.add('/home/data/some/other/path')
        then:
        tree.longest().size() ==1
        tree.longest().get(0) == '/home/data'

    }

    def testOnlyRoot() {

        when:
        def tree = new PathTrie()
        tree.add('/')
        then:
        tree.longest().size() == 0

        when:
        tree = new PathTrie()
        tree.add('a')
        then:
        tree.longest().size() == 1
        tree.longest().get(0) == 'a'

        when:
        tree = new PathTrie()
        tree.add('/a')
        then:
        tree.longest().size() == 1
        tree.longest().get(0) == '/a'
    }


}
