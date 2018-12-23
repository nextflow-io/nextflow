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

