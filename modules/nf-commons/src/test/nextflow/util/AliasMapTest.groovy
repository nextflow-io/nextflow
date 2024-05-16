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
 */

package nextflow.util

import spock.lang.Specification
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class AliasMapTest extends Specification {

    def 'should convert keys from kebab case to camel case' () {

        when:
        def map = new AliasMap()
        map['alphaBeta'] = 1
        map['alpha-beta'] = 10
        then:
        map['alphaBeta'] == 10
        !map.containsKey('alpha-beta')

        when:
        map = new AliasMap()
        map['aaa-bbb-ccc'] = 1
        map['aaaBbbCcc'] = 10
        then:
        map['aaaBbbCcc'] == 10
        !map.containsKey('aaa-bbb-ccc')

    }

    def 'should convert kebab case to camel case' () {

        expect:
        AliasMap.kebabToCamelCase('a') == 'a'
        AliasMap.kebabToCamelCase('A') == 'A'
        AliasMap.kebabToCamelCase('a-b-c-') == 'aBC'
        AliasMap.kebabToCamelCase('aa-bb-cc') == 'aaBbCc'
        AliasMap.kebabToCamelCase('alpha-beta-delta') == 'alphaBetaDelta'
        AliasMap.kebabToCamelCase('Alpha-Beta-delta') == 'AlphaBetaDelta'
    }

    def 'should convert camel case to kebab case' () {

        expect:
        AliasMap.camelToKebabCase('alphaBetaDelta') == 'alpha-beta-delta'
        AliasMap.camelToKebabCase('AlphaBetaDelta') == 'Alpha-beta-delta'
        AliasMap.camelToKebabCase('Field1') == 'Field1'
        AliasMap.camelToKebabCase('FieldUno') == 'Field-uno'
        AliasMap.camelToKebabCase('FieldUNO') == 'Field-UNO'
        AliasMap.camelToKebabCase('FieldA') == 'Field-A'
        AliasMap.camelToKebabCase('FieldAB') == 'Field-AB'
        AliasMap.camelToKebabCase('FieldAb') == 'Field-ab'
    }

}