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

package nextflow.ui


import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TextLabelTest extends Specification {

    def "test pad" () {

        expect:
        TextLabel.of('Hello').left().width(8).toString() == "Hello   "
        TextLabel.of('Hello').right().width(8).toString() == "   Hello"
        TextLabel.of('1234').number().width(8).toString() == "    1234"

    }

    def 'text align' () {

        /*
         * numbers are aligned automatically to the right
         */
        expect:
        new TextLabel('Hello').width(8).toString() == "Hello   "
        new TextLabel('1234').width(8).toString()  == "    1234"
        new TextLabel(1234).width(8).toString()    == "    1234"
        new TextLabel(1234).width(8).left().toString() == "1234    "

    }

    def 'text max' () {
        expect:
        new TextLabel('Hello world!   ').max(12).toString() == 'Hello world!'
        new TextLabel('Hello world!').max(12).toString() == 'Hello world!'
        new TextLabel('Hello world!').max(10).toString() == 'Hello wo..'
        new TextLabel('Hello world!').max(8).toString()  == 'Hello ..'
        new TextLabel('Hello world!').max(8).width(20).toString()  == 'Hello ..'

        new TextLabel(12345).max(6).toString() == '12345'
        new TextLabel(12345).max(5).toString() == '12345'
        new TextLabel(12345).max(4).toString() == '..45'

    }


    def "test apply decorator " () {

        when:
        def label = TextLabel.of( 'The cat eat the food' )
        def deco  = { obj, val -> "* ${val} *".toString() } as LabelDecorator
        label << deco

        then:
        label.toString() == "The cat eat the food"
        label.switchOn().toString() == "* The cat eat the food *"

    }

//    def "test add decorator" () {
//
//        when:
//        def label = TextLabel.of('The Fact cat') << AnsiStyle.style().bold()
//
//        then:
//        label.decorators.find { it instanceof AnsiStyle }
//        label.decorators.size() == 1
//
//    }

}