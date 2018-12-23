/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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
