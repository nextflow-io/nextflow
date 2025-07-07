/*
 * Copyright 2013-2024, Seqera Labs
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
        node.addNode('db')
        node.addNode('data')
        then:
        node.getChildren().size() == 2
        node.getChildren() *. node == ['db', 'data']

    }

    def testAppendTwice() {
        when:
        def node = new Trie<String>('/')
        node.addNode('db')
        node.addNode('db')
        then:
        node.getChildren().size() == 1
        node.getChildren() *. node == ['db']
    }

    def testAppendConcat() {
        when:
        def node = new Trie<String>('/')
        node.addNode('db').addNode('some').addNode('path')
        then:
        node.children.size() == 1
        node.children[0].node == 'db'
        node.children[0].children[0].node == 'some'
        node.children[0].children[0].children[0].node == 'path'
        node.children[0].children[0].children[0].children == null
    }


    def 'should add path'() {
        when:
        def node = new Trie<String>('/')
        node.addPath( 'alpha', 'beta' )
        then:
        node.children.size() == 1
        and:
        node.children[0].node == 'alpha'
        node.children[0].children[0].node == 'beta'
        node.children[0].children[0].children == null

        when:
        node.addPath('alpha','beta','delta','omega')
        then:
        node.children[0].children[0].node == 'beta'
        node.children[0].children[0].children[0].node == 'delta'
        node.children[0].children[0].children[0].children[0].node == 'omega'

    }

    def testAppend2() {
        when:
        def node = new Trie<String>('/')
        node.addPath( 'db', 'some', 'file.txt' )
        node.addPath( 'db', 'some', 'file.fa' )

        then:
        node.children.size() == 1
        node.children[0].node == 'db'

        node.children[0].children.size() == 1
        node.children[0].children[0].node == 'some'

        node.children[0].children[0].children.size() == 2
        node.children[0].children[0].children[0].node == 'file.txt'
        node.children[0].children[0].children[0].children == null

        node.children[0].children[0].children[1].node == 'file.fa'
        node.children[0].children[0].children[1].children == null
    }

    def testLongest() {

        when:
        def node = new Trie<String>('db')
        then:
        node.longest() == ['db']

        when:
        node = new Trie<String>('db')
        node.addPath('data','path','file1')
        node.addPath('data','path','file2')
        then:
        node.longest() == ['db','data','path']

    }

    def 'should traverse tree' () {
        given:
        def trie = new Trie<String>('root')
        and:
        trie.addPath('alpha','beta','delta')
        trie.addPath('alpha','beta','gamma')
        trie.addPath('alpha','beta','gamma','omega')
        trie.addPath('alpha','episol','teta')
        when:
        def result = trie.traverse()
        then:
        result.size() == 3
        and:
        result[0].join(', ') == 'root, alpha, beta, delta'
        result[1].join(', ') == 'root, alpha, beta, gamma, omega'
        result[2].join(', ') == 'root, alpha, episol, teta'
    }

}

