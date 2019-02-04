/*
 * Copyright 2013-2019, Centre for Genomic Regulation (CRG)
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
