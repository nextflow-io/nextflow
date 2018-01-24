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
class TrieTest extends Specification {

    def testTrie() {

        when:
        def node = new Trie<String>('/')
        node.append('db')
        node.append('data')
        then:
        node.getChildren().size() == 2
        node.getChildren() *. vertex == ['db','data']

    }

    def testAppendTwice() {
        when:
        def node = new Trie<String>('/')
        node.append('db')
        node.append('db')
        then:
        node.getChildren().size() == 1
        node.getChildren() *. vertex == ['db']
    }

    def testAppendConcat() {
        when:
        def node = new Trie<String>('/')
        node.append('db').append('some').append('path')
        then:
        node.children.size() == 1
        node.children[0].vertex == 'db'
        node.children[0].children[0].vertex == 'some'
        node.children[0].children[0].children[0].vertex == 'path'
        node.children[0].children[0].children[0].children == null
    }


    def testAppendList() {
        when:
        def node = new Trie<String>('/')
        node.append( 'db', 'some', 'path' )
        then:
        node.children.size() == 1
        node.children[0].vertex == 'db'
        node.children[0].children[0].vertex == 'some'
        node.children[0].children[0].children[0].vertex == 'path'
        node.children[0].children[0].children[0].children == null
    }

    def testAppend2() {
        when:
        def node = new Trie<String>('/')
        node.append( 'db', 'some', 'file.txt' )
        node.append( 'db', 'some', 'file.fa' )

        then:
        node.children.size() == 1
        node.children[0].vertex == 'db'

        node.children[0].children.size() == 1
        node.children[0].children[0].vertex == 'some'

        node.children[0].children[0].children.size() == 2
        node.children[0].children[0].children[0].vertex == 'file.txt'
        node.children[0].children[0].children[0].children == null

        node.children[0].children[0].children[1].vertex == 'file.fa'
        node.children[0].children[0].children[1].children == null
    }

    def testLongest() {

        when:
        def node = new Trie<String>('db')
        then:
        node.longest() == ['db']

        when:
        node = new Trie<String>('db')
        node.append('data','path','file1')
        node.append('data','path','file2')
        then:
        node.longest() == ['db','data','path']

    }

}

