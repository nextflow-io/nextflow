/*
 * Copyright 2013-2023, Seqera Labs
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
 *
 */

package io.seqera.wave.util

import spock.lang.Specification

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class TemplateRendererTest extends Specification {

    def 'should replace vars' () {
        given:
        def binding = [foo: 'Hello', bar: 'world']

        expect:
        TemplateRenderer.replace0('{{foo}}', binding) == 'Hello'
        TemplateRenderer.replace0('{{foo}} ', binding) == 'Hello '
        TemplateRenderer.replace0('{{foo}}\n', binding) == 'Hello\n'
        TemplateRenderer.replace0('  {{foo}}', binding) == '  Hello'
        TemplateRenderer.replace0('  {{foo}}\n', binding) == '  Hello\n'
        TemplateRenderer.replace0('  ${foo}', binding)  == '  ${foo}'
        TemplateRenderer.replace0('  ${{foo}}', binding) == '  ${{foo}}'
        TemplateRenderer.replace0('{{foo}}', [foo:'']) == ''
        TemplateRenderer.replace0('{{foo}}', [foo:null]) == null
        TemplateRenderer.replace0('  {{foo}}\n', [foo:null]) == null
        TemplateRenderer.replace0('', binding) == ''
        TemplateRenderer.replace0(null, binding) == null

        TemplateRenderer.replace0('{{foo}} {{bar}}!', binding) == 'Hello world!'
        TemplateRenderer.replace0('abc {{foo}} pq {{bar}} xyz', binding) == 'abc Hello pq world xyz'
        TemplateRenderer.replace0('{{foo}} 123 {{bar}} xyz {{foo}}', binding) == 'Hello 123 world xyz Hello'
        TemplateRenderer.replace0('1{{foo}}2{{foo}}3', [foo:'']) == '123'
        TemplateRenderer.replace0('1{{foo}}2{{foo}}3', [foo:null]) == '123'

        when:
        TemplateRenderer.replace0('{{x1}}', binding)
        then:
        def e = thrown(IllegalArgumentException)
        e.message == 'Missing template key: x1'

        when:
        TemplateRenderer.replace0('{{foo}} {{x2}}', binding)
        then:
        e = thrown(IllegalArgumentException)
        e.message == 'Missing template key: x2'

    }

    def 'should render template' () {
        given:
        def template = "Hello, {{name}}!\n" +
                "Today is {{day}} and the weather is {{weather}}.";
        and:
        def binding = new HashMap<String,String>();
        binding.put("name", "John");
        binding.put("day", "Monday");
        binding.put("weather", "sunny");

        when:
        def result = TemplateRenderer.render(template, binding);

        then:
        result == 'Hello, John!\nToday is Monday and the weather is sunny.'
    }

    def 'should render a template with comment'() {
        given:
        def template = """\
            ## remove this comment
            1: {{alpha}}
            2: {{delta}} {{delta}}
            3: {{gamma}} {{gamma}} {{gamma}} 
            4: end
            """.stripIndent()
        and:
        def binding = new HashMap<String,String>();
        binding.put("alpha", "one");
        binding.put("delta", "two");
        binding.put("gamma", "three");

        when:
        def result = TemplateRenderer.render(new ByteArrayInputStream(template.bytes), binding);

        then:
        result == """\
            1: one
            2: two two
            3: three three three 
            4: end
            """.stripIndent()
    }


    def 'should render a template using an input stream'() {
        given:
        def template = """\
            {{one}}
            {{two}}
            xxx
            {{three}}            
            zzz
            """.stripIndent()
        and:
        def binding = [
                one: '1',   // this is rendered
                two:null,   // a line containing a null variable is not rendered
                three:''    // empty value is considered ok
        ]

        when:
        def result = TemplateRenderer.render(new ByteArrayInputStream(template.bytes), binding);

        then:
        result == """\
            1
            xxx
            
            zzz
            """.stripIndent()
    }

    def 'should render template with indentations' () {
        given:
        def binding = [foo: 'Hello', bar: 'world']

        when:
        def result = TemplateRenderer.render('{{foo}}\n{{bar}}', binding)
        then:
        result == 'Hello\nworld'

        when:
        def template = '''\
              {{foo}}
            {{bar}}
            '''.stripIndent()
        result = TemplateRenderer.render(template, [foo:'11\n22\n33', bar:'Hello world'])
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
        result = TemplateRenderer.render(template, [x1:'aa\nbb\n', x2:null, x3:'pp\nqq'])
        then:
        result == '''\
            aa
            bb
            
            pp
            qq
            '''.stripIndent()

    }

}
