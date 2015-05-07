/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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
class EscapeTemplateEngineTest extends Specification {

    def 'should not interpolate dollar prefixed variables'() {

        given:
        def text = '''
        alpha: $bash_var
        delta: $(bash expr)
        gamma: ${bash var}
        groovy1: &foo
        groovy2: &obj.beta
        groovy3: &{obj.pico}
        groovy4: &{(x + y) / z}
        groovy5: '&{(x + z) / y}'
        groovy6: "&{x + y / z}"
        do not resolve escape: \t \n
        ignore this: &!
        .. and this: & hola
        .. and this: &\t
        '''
        .stripIndent()
//        theta: &! ignore this

        when:
        def binding = [foo: 'bar', x: 9, y: 1, z: 2, obj: [beta: 'Hello', pico: 'World']]
        def template = new EscapeTemplateEngine().createTemplate(text)
        then:
        template.make(binding).toString() == '''
            alpha: $bash_var
            delta: $(bash expr)
            gamma: ${bash var}
            groovy1: bar
            groovy2: Hello
            groovy3: World
            groovy4: 5
            groovy5: '11'
            groovy6: "9.5"
            do not resolve escape: \t \n
            ignore this: &!
            .. and this: & hola
            .. and this: &\t
            '''
            .stripIndent()
    }


    def 'should interpolate file template' () {

        given:
        def template = this.class.getResourceAsStream('template.txt')
        when:
        def binding = [foo: 'bar', x: 9, y: 1, z: 2, obj: [beta: 'Hello', pico: 'World']]
        def result = new EscapeTemplateEngine().createTemplate(new InputStreamReader(template)).make(binding)
        then:
        result.toString() == new InputStreamReader(this.class.getResourceAsStream('template.result')).text

    }
}
