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
class MustacheTemplateEngineTest extends Specification {

    def 'should match variable placeholder' () {

        when:
        def m1 = MustacheTemplateEngine.VAR1.matcher('{{FOO}}')
        then:
        m1.matches()
        m1.groupCount() == 3
        m1.group(1) == ''
        m1.group(2) == 'FOO'
        m1.group( 3) == ''

        when:
        def m2 = MustacheTemplateEngine.VAR1.matcher('  {{this_123}} ')
        then:
        m2.matches()
        m2.groupCount() == 3
        m2.group(1) == '  '
        m2.group(2) == 'this_123'
        m2.group( 3) == ' '

        when:
        def m3 = MustacheTemplateEngine.VAR1.matcher('{{that-abc}}')
        then:
        m3.matches()
        m3.groupCount() == 3
        m3.group(1) == ''
        m3.group(2) == 'that-abc'
        m3.group( 3) == ''

        expect:
        !MustacheTemplateEngine.VAR1.matcher('{{}}').matches()
        !MustacheTemplateEngine.VAR1.matcher('${{FOO}}').matches()
        !MustacheTemplateEngine.VAR1.matcher('${FOO}').matches()
        !MustacheTemplateEngine.VAR1.matcher(' ${FOO}').matches()
        !MustacheTemplateEngine.VAR1.matcher('xxx{{FOO}}').matches()
    }

    def 'should find var placeholders' () {

        when:
        def m1 = MustacheTemplateEngine.VAR2.matcher('{{FOO}}')
        then:
        m1.find()
        m1.groupCount() == 1
        [ m1.start(1), m1.end(1) ] == [2,5]


        when:
        def m2 = MustacheTemplateEngine.VAR2.matcher('abc {{FOO}} {{BAR}} xyz')
        then:
        m2.find()
        m2.groupCount() == 1
        [m2.start(1), m2.end(1)] == [6,9]

        when:
        def m3 = MustacheTemplateEngine.VAR2.matcher('${{FOO}}')
        then:
        !m3.find()
    }

    def 'should replace vars' () {
        given:
        def binding = [foo: 'Hello', bar: 'world']
        def mustache = new MustacheTemplateEngine()

        expect:
        mustache.replace0('{{foo}}', binding) == 'Hello'
        mustache.replace0('  {{foo}}', binding) == '  Hello'
        mustache.replace0('  ${foo}', binding)  == '  ${foo}'
        mustache.replace0('  ${{foo}}', binding) == '  ${{foo}}'
        mustache.replace0('{{foo}}', [foo:'']) == null
        mustache.replace0('  {{foo}}', [foo:null]) == null
        mustache.replace0('', binding) == ''
        mustache.replace0(null, binding) == null

        mustache.replace0('{{foo}} {{bar}}!', binding) == 'Hello world!'
        mustache.replace0('abc {{foo}} pq {{bar}} xyz', binding) == 'abc Hello pq world xyz'
        mustache.replace0('{{foo}} 123 {{bar}} xyz {{foo}}', binding) == 'Hello 123 world xyz Hello'
        mustache.replace0('1{{foo}}2{{foo}}3', [foo:'']) == '123'

        when:
        mustache.replace0('{{x1}}', binding)
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Missing template key: x1'

        when:
        mustache.replace0('{{foo}} {{x2}}', binding)
        then:
        e = thrown(IllegalArgumentException)
        e.message == 'Missing template key: x2'

    }

    def 'should render template' () {

        given:
        def binding = [foo: 'Hello', bar: 'world']
        def mustache = new MustacheTemplateEngine()

        when:
        def result = mustache.render('{{foo}}\n{{bar}}', binding)
        then:
        result == 'Hello\nworld\n'

        when:
        def template = '''\
              {{foo}}
            {{bar}}
            '''.stripIndent()
        result = mustache.render(template, [foo:'11\n22\n33', bar:'Hello world'])
        then:
        result == '''\
              11
              22
              33
            Hello world
            '''.stripIndent()


        when:
        template = '''\
            {{x1}}
            {{x2}}
            {{x3}}
            '''.stripIndent()
        result = mustache.render(template, [x1:'aa\nbb\n', x2:null, x3:'pp\nqq'])
        then:
        result == '''\
            aa
            bb
            pp
            qq
            '''.stripIndent()

    }

    def 'should render empty lines' () {
        given:
        def binding = [foo: 'Hello', bar: 'world']
        def mustache = new MustacheTemplateEngine()

        when:
        def template = '''\
            {{foo}}
            
            {{bar}}!
            '''.stripIndent()
        def result = mustache.render(template, binding)
        then:
        result == '''\
                Hello
                
                world!
                '''.stripIndent()
    }

    def 'should strip comments from template' () {

        given:
        def engine = new MustacheTemplateEngine() {
            @Override protected boolean accept(String line) { line.endsWith('2')}
        }
        def template = '''\
            line 1
            line 2
            line 3
            '''.stripIndent()

        expect:
        engine.render(template, [:]) == '''\
                line 2
                '''.stripIndent()

    }


}
