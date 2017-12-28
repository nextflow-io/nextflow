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

package nextflow.splitter

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CsvSplitterTest extends Specification {

    def text = '''
        alpha,beta,delta
        gamma,,zeta
        eta,theta,iota
        mu,nu,xi
        pi,rho,sigma
        '''
        .stripIndent().trim()

    def testSplitRows() {

        when:
        def items = new CsvSplitter().target(text).list()

        then:
        items.size() == 5
        items[0] instanceof List
        items[0] == ['alpha', 'beta', 'delta']
        items[1] == ['gamma', '', 'zeta']
        items[2] == ['eta', 'theta', 'iota']
        items[3] == ['mu', 'nu', 'xi']
        items[4] == ['pi', 'rho', 'sigma']
    }

    def testSplitRowsWithLimit() {

        when:
        def items = new CsvSplitter(limit: 3).target(text).list()

        then:
        items.size() == 3
        items[0] instanceof List
        items[0] == ['alpha', 'beta', 'delta']
        items[1] == ['gamma', '', 'zeta']
        items[2] == ['eta', 'theta', 'iota']

    }


    def testSkipRows() {
        when:
        def items = new CsvSplitter().target(text).options(skip:3).list()
        then:
        items.size() == 2
        items[0] instanceof List
        items[0] == ['mu', 'nu', 'xi']
        items[1] == ['pi', 'rho', 'sigma']

        when:
        items = new CsvSplitter().target(text).options(skip:1, limit: 3).list()
        then:
        items.size() == 3
        items[0] instanceof List
        items[0] == ['gamma', '', 'zeta']
        items[1] == ['eta', 'theta', 'iota']
        items[2] == ['mu', 'nu', 'xi']
    }

    def testSplitWithCount() {

        when:
        def groups = new CsvSplitter(by:3).target(text).list()

        then:
        groups.size() == 2

        groups[0][0] == ['alpha', 'beta', 'delta']
        groups[0][1] == ['gamma', '', 'zeta']
        groups[0][2] == ['eta', 'theta', 'iota']

        groups[1][0] == ['mu', 'nu', 'xi']
        groups[1][1] == ['pi', 'rho', 'sigma']

    }

    def testStrip( ) {

        when:
        def splitter = new CsvSplitter().options(quote: '"')
        then:
        splitter.strip('hello') == 'hello'
        splitter.strip('"hello"') == 'hello'
        splitter.strip('\'hello\'') == '\'hello\''
        splitter.strip('""') == ''
        splitter.strip('"') == '"'

        when:
        splitter = new CsvSplitter().options(quote: '"', strip: true)
        then:
        splitter.stripBlanks
        splitter.strip('hello ') == 'hello'
        splitter.strip('"hello"') == 'hello'
        splitter.strip('" hello "  ') == ' hello '
        splitter.strip(' " " ') == ' '
        splitter.strip(' "" ') == ''

    }

    def testSplitCsvWithHeaderCols() {

        when:
        def items = new CsvSplitter().target(text).options(header:['col1','col2','col3']).list()
        then:
        items.size() == 5

        items[0].col1 == 'alpha'
        items[0].col2 == 'beta'
        items[0].col3 == 'delta'

        items[1].col1 == 'gamma'
        items[1].col2 == ''
        items[1].col3 == 'zeta'

        items[2].col1 == 'eta'
        items[2].col2 == 'theta'
        items[2].col3 == 'iota'
    }

    def testSplitCsvWithHeaderLine() {

        when:
        def items = new CsvSplitter().target('x,y,z\n' + text).options(header:true).list()
        then:
        items.size() == 5

        items[0].x == 'alpha'
        items[0].y == 'beta'
        items[0].z == 'delta'

        items[1].x == 'gamma'
        items[1].y == ''
        items[1].z == 'zeta'

        items[2].x == 'eta'
        items[2].y == 'theta'
        items[2].z == 'iota'
    }

    def testSplitCsvGroupMap() {

        when:
        def items = new CsvSplitter().target(text).options(by:3, header:['x','y','z']).list()
        then:
        items.size() == 2
        items[0] == [ [x:'alpha',y:'beta',z:'delta'], [x:'gamma',y:'',z:'zeta'], [x:'eta',y:'theta',z:'iota'] ]
        items[1] == [ [x:'mu',y:'nu',z:'xi'], [x:'pi', y:'rho', z:'sigma'] ]

    }

    def testIllegalRecordMode() {

        when:
        new CsvSplitter().options(record:true)
        then:
        thrown(IllegalArgumentException)

    }

    def testIllegalFileMode() {

        when:
        new CsvSplitter().options(file:true)
        then:
        thrown(IllegalArgumentException)

    }

}
