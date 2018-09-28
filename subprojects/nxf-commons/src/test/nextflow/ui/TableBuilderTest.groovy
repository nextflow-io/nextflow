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

package nextflow.ui

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TableBuilderTest extends Specification {

    def 'test trackSize' () {

        setup:
        def table = new TableBuilder()

        when:
        table.append([1,22,333])
        def row1 = new ArrayList(table.dim)

        table.append(['abc','x','y','12345', 'z'])
        def row2 = new ArrayList(table.dim)

        then:
        row1 == [1,2,3]
        row2 == [3,2,3,5,1]

    }

    def 'test table builder' () {

        when:
        def table = new TableBuilder()
        table.append([ 1, 2, 3 ] )
        table.append([ 4, 5, 6 ] )

        then:
        table.toString() == "1 2 3\n4 5 6"

    }



    def 'test table builder with header' () {

        when:
        def table = new TableBuilder().setHeaders('alpha', 'beta', 'delta')
        table.append([ 1, 2, 3 ] )
        table.append([ 4, 5, 6 ] )

        then:
        table.toString() == """
        alpha beta delta
            1    2     3
            4    5     6
        """
                .stripIndent().trim()

    }


    def 'test table build with long string' () {

        when:
        def table = new TableBuilder().setHeaders('alpha', 'beta', 'delta__')
        table.append( 'hola', '1', 'some long string' )
        table.append( 'ciao', '2', 'short string' )
        table.append( 'hello', '3', 'more text' )
        table.setMaxColsWidth(5,5,7)


        then:
        table.toString() == """
        alpha beta delta__
        hola     1 some ..
        ciao     2 short..
        hello    3 more ..
        """
                .stripIndent().trim()

    }

    def 'test table build with shit ' () {

        when:
        def table = new TableBuilder()
        table.head('alpha').head('beta').head('delta')

        table << 1 << 2 << 3 << table.closeRow()
        table << 4 << 5 << 6

        then:
        table.toString() == """
        alpha beta delta
            1    2     3
            4    5     6
        """
                .stripIndent().trim()

    }

    def 'test table head with max' () {

        when:
        def table = new TableBuilder()
        table.head('x',3).head('y____',5)

        table << 'a' << '11111111' << table.closeRow()
        table << 'bb' << '22222222' << table.closeRow()
        table << 'ccc' << '33333333' << table.closeRow()

        then:
        table.toString() == """
        x   y____
        a   ..111
        bb  ..222
        ccc ..333
        """
                .stripIndent().trim()

    }

    def 'test table head with align' () {

        when:
        def table = new TableBuilder()
                .head('x__',TextLabel.Align.RIGHT)
                .head('y__', TextLabel.Align.LEFT)

        table << 'a'  << '1' << table.closeRow()
        table << 'bb' << '22' << table.closeRow()
        table << '9'  << '333' << table.closeRow()

        then:
        table.toString() == """
        x__ y__
          a 1  \n\
         bb 22 \n\
          9 333
        """
                .stripIndent().trim()

    }

}