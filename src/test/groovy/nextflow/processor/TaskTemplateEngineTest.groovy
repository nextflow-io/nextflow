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

package nextflow.processor

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TaskTemplateEngineTest extends Specification {

    def 'should parse GString like template' () {

        expect:
        new TaskTemplateEngine().render('Hello') == 'Hello'
        new TaskTemplateEngine().render('$a ${b}', [a: 'Hello', b: 'world!']).toString() == 'Hello world!'
        new TaskTemplateEngine().render('$obj.a \n ${obj.b}', [obj: [a: 'Hello', b: 'world!']]).toString() == 'Hello \n world!'

    }

    def 'should ignore dollar variable' () {

        given:
        def text = '''
        alpha: $bash_var
        delta: $(bash expr)
        gamma: ${bash var}
        say..: !{foo} !{obj.bar}
        not..: !Hello
        $XXXXXX!{ZZZ}
        '''
                .stripIndent()

        when:
        def binding = [foo: 'Hello', obj: [bar: 'world!'], ZZZ: '___']
        def result = new TaskTemplateEngine()
                    .setPlaceholder('!' as char)
                    .createTemplate(text)
                    .make(binding)
                    .toString()
        then:
        result == '''
            alpha: $bash_var
            delta: $(bash expr)
            gamma: ${bash var}
            say..: Hello world!
            not..: !Hello
            $XXXXXX___
            '''
                .stripIndent()

    }

    def 'should not interpolate dollar prefixed variables'() {

        given:
        def text = '''
        alpha: $bash_var
        delta: $(bash expr)
        gamma: ${bash var}
        groovy1: !foo
        groovy2: !obj.beta
        groovy3: !{obj.pico}
        groovy4: !{(x + y) / z}
        groovy5: '!{(x + z) / y}'
        groovy6: "!{x + y / z}"
        do not resolve escape: \t \n
        ignore this: !&
        .. and this: ! hola
        .. and this: !\t
        '''
        .stripIndent()

        when:
        def binding = [foo: 'bar', x: 9, y: 1, z: 2, obj: [beta: 'Hello', pico: 'World']]
        def result = new TaskTemplateEngine(placeholder: '!' as char, enableShortNotation: true) .render(text, binding)
        then:
        result == '''
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
            ignore this: !&
            .. and this: ! hola
            .. and this: !\t
            '''
            .stripIndent()
    }

    def 'should use custom placeholder' () {

        given:
        def text = '''
        alpha: $bash_var
        delta: $(bash expr)
        gamma: ${bash var}
        say..: #{foo} #{obj.bar}
        not..: #Hello
        '''
                .stripIndent()

        when:
        def binding = [foo: 'Hello', obj: [bar: 'world!']]
        def result = new TaskTemplateEngine(placeholder: '#' as char).render(text, binding)
        then:
        result == '''
            alpha: $bash_var
            delta: $(bash expr)
            gamma: ${bash var}
            say..: Hello world!
            not..: #Hello
            '''
                .stripIndent()
    }

    def 'should interpolate file template' () {

        given:
        def template = this.class.getResourceAsStream('/nextflow/processor/test-file.tpl')
        when:
        def binding = [foo: 'bar', x: 9, y: 1, z: 2, obj: [beta: 'Hello', pico: 'World']]
        def result = new TaskTemplateEngine()
                    .setPlaceholder('!' as char)
                    .createTemplate(new InputStreamReader(template))
                    .make(binding)
                    .toString()
        then:
        result == new InputStreamReader(this.class.getResourceAsStream('/nextflow/processor/test-result.tpl')).text

    }

    def 'should not change the binding map object content' () {

        given:
        def binding = [x: 'Hello', y: 0]
        when:
        def original = new HashMap(binding)
        def result = new TaskTemplateEngine().render('Say: $x ${++y} and ${++y} times', binding)
        then:
        result == 'Say: Hello 1 and 2 times'
        original.size() == binding.size()

    }

    def 'should hold a print writer object' () {

        given:
        def writer = new StringWriter()
        def context = new TaskTemplateEngine.TemplateBinding([x:'foo', z:'bar'], writer)

        when:
        context.getVariable('__$$_out') << 'Hello world'

        then:
        writer.toString() == 'Hello world'
        context.getVariable('x') == 'foo'
        context.getVariable('z') == 'bar'
        context.x == 'foo'
        context.z == 'bar'
        context.hasVariable('x')
        context.hasVariable('z')
        !context.hasVariable('y')

    }
}
