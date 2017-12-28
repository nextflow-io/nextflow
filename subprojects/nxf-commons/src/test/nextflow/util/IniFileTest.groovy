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
class IniFileTest extends Specification {

    def 'test read ini file' () {

        given:
        def text = '''
            [eng]
            x = Hello
            y = world
            z = 123
            b = true

            [esp]
            x = Hola
            y = mundo
            z = 1.23

            [ita]
            x = Ciao
            y = mondo
            z = 3.1415926535
            '''

        when:
        def ini = new IniFile().load(new StringReader(text))
        then:
        ini.getString('eng','x') == 'Hello'
        ini.getString('eng','y') == 'world'
        ini.getInt('eng','z') == 123
        ini.getBool('eng', 'b') == true

        ini.getString('esp','x') == 'Hola'
        ini.getString('esp','y') == 'mundo'
        ini.getDouble('esp','z') == 1.23d

        ini.getString('ita','x') == 'Ciao'
        ini.getString('ita','y') == 'mondo'
        ini.getFloat('ita','z') == 3.1415926535f

    }

    def 'default values' () {
        given:
        def text = '''
            [eng]
            x = Hello
            y = world
            '''
        when:
        def ini = new IniFile().load(new StringReader(text))
        then:
        ini.getString('eng', 'q') == null
        ini.getString('eng', 'q', 'Missing !!') == 'Missing !!'
        ini.getString('xxx', 'q', 'Missing !!') == 'Missing !!'

        ini.getInt('eng', 'q') == 0
        ini.getInt('eng', 'q', -1) == -1
        ini.getInt('xxx', 'q', -1) == -1

        ini.getFloat('eng', 'q') == 0
        ini.getFloat('eng', 'q', -2) == -2
        ini.getFloat('xxx', 'q', -2) == -2

        ini.getDouble('eng', 'q') == 0
        ini.getDouble('eng', 'q', -3) == -3
        ini.getDouble('xxx', 'q', -3) == -3

        ini.getBool('eng', 'q') == false
        ini.getBool('eng', 'q', true) == true
        ini.getBool('xxx', 'q', true) == true
    }

    def 'test dynamic' () {

        given:
        def text = '''
            [eng]
            x = Hello
            y = world
            z = 123

            [esp]
            x = Hola
            y = mundo
            z = 1.23

            [ita]
            x = Ciao
            y = mondo
            z = 3.1415926535
            '''

        when:
        def ini = new IniFile().load(new StringReader(text))

        then:
        ini.eng.x == 'Hello'
        ini.eng.y == 'world'

    }

}
