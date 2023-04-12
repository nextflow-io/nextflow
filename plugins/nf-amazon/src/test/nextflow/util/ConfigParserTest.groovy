/*
 * Copyright 2020-2022, Seqera Labs
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

package nextflow.util

import nextflow.cloud.aws.util.ConfigParser
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class ConfigParserTest extends Specification {

    def 'should parse section' () {
        given:
        def parser = new ConfigParser()

        expect:
        parser.parseSection(LINE) == EXPECTED
        
        where:
        LINE            | EXPECTED
        'foo'           | null
        '[foo'          | null
        and:
        '[foo]'         | 'foo'
        '[profile foo]' | 'foo'
    }

    def 'should parse config' () {
        given:
        def parser = new ConfigParser()
        def CONFIG = '''
        [foo]
        one = 1 
        two = 2 
        [bar]
        alpha = 3
        gamma = 4        
        '''.stripIndent()

        when:
        parser.parseConfig(CONFIG)

        then:
        parser.content.size() == 2
        parser.content['foo']  == ['one = 1', 'two = 2']
        parser.content['bar']  == ['alpha = 3', 'gamma = 4']
    }

    def 'should not merge overlapping keys' () {
        given:
        def parser = new ConfigParser()
        def CONFIG1 = '''
        [alpha]
        a1=1
        [beta]
        b2=2
        b3=3
        '''.stripIndent()

        def CONFIG2 = '''
        [beta]
        b3=30
        b4=4
        '''.stripIndent()


        when:
        parser.parseConfig(CONFIG1)
        parser.parseConfig(CONFIG2)

        then:
        parser.content.size() == 2
        and:
        parser.content['alpha'] == ['a1=1']
        parser.content['beta'] ==  ['b2=2','b3=3','b4=4']

    }

    def 'should load and merge config' () {
        given:
        def parser = new ConfigParser()
        def CONFIG1 = '''
        [alpha]
        a1
        '''.stripIndent()

        def CONFIG2 = '''
        [beta]
        b1
        '''.stripIndent()

        def CONFIG3 = '''
        [alpha]
        a2
        
        [beta]
        b2
        
        [omega]
        z9        
        '''.stripIndent()

        when:
        parser.parseConfig(CONFIG1)
        parser.parseConfig(CONFIG2)
        parser.parseConfig(CONFIG3)

        then:
        parser.content.size() == 3
        and:
        parser.content['alpha'] == ['a1','a2']
        parser.content['beta'] ==  ['b1','b2']
        parser.content['omega'] ==  ['z9']

        expect:
        parser.text() == '''\
            [alpha]
            a1
            a2
            [beta]
            b1
            b2
            [omega]
            z9
            '''.stripIndent()

    }


    @Unroll
    def 'should match key' () {
        given:
        def parser = new ConfigParser()

        expect:
        parser.findKey(LINE) == EXPECTED

        where:
        LINE        | EXPECTED
        'foo'       | null
        'foo='      | 'foo'
        'foo=1'     | 'foo'
        ' foo = 1 ' | 'foo'
        ' foo  =1 ' | 'foo'

    }
}
