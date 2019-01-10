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
