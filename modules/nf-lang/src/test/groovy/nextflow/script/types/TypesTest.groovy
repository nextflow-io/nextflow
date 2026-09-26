/*
 * Copyright 2024-2025, Seqera Labs
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

package nextflow.script.types

import nextflow.script.dsl.Types
import org.codehaus.groovy.ast.ClassHelper
import org.codehaus.groovy.ast.ClassNode
import org.codehaus.groovy.ast.GenericsType
import spock.lang.Specification

/**
 *
 * @author Ben Sherman <bentshermann@gmail.com>
 */
class TypesTest extends Specification {

    def 'should render a type' () {
        when:
        def cn = ClassHelper.MAP_TYPE
        then:
        Types.getName(cn) == 'Map<K, V>'

        when:
        cn = new ClassNode(Map.class)
        then:
        Types.getName(cn) == 'Map'

        when:
        cn.setGenericsTypes(new GenericsType[] {
            new GenericsType(ClassHelper.STRING_TYPE),
            new GenericsType(ClassHelper.STRING_TYPE)
        })
        then:
        Types.getName(cn) == 'Map<String, String>'
    }

    static class Sample implements Record {}

    def 'should render a runtime type' () {
        expect:
        Types.getName(TYPE) == NAME

        where:
        TYPE                                   | NAME
        java.nio.file.Paths.get('/a').class    | 'Path'
        ArrayList                              | 'List'
        Sample                                 | 'Sample'
        Record                                 | 'Record'
        String                                 | 'String'
    }

    def 'should render the return type of a method' () {
        when:
        def cn = ClassHelper.makeCached(TaskConfig)
        def returnType = cn.getDeclaredMethods('getExt')[0].getReturnType()
        then:
        Types.getName(returnType) == 'Map<String, String>'
    }

    def 'should render a functional type' () {
        given:
        def cn = ClassHelper.makeCached(Channel)

        when:
        def param = cn.getDeclaredMethods('filter')[0].getParameters()[0]
        then:
        Types.getName(param.getType()) == '(E) -> Boolean'

        when:
        param = cn.getDeclaredMethods('map')[0].getParameters()[0]
        then:
        Types.getName(param.getType()) == '(E) -> R'

        when:
        param = cn.getDeclaredMethods('reduce')[0].getParameters().last()
        then:
        Types.getName(param.getType()) in ['(E, E) -> E', '(R, E) -> R']
    }

    def 'should determine whether a type is assignable to another type' () {
        expect:
        Types.isAssignableFrom(ClassHelper.makeCached(TARGET), ClassHelper.makeCached(SOURCE)) == RESULT

        where:
        TARGET      | SOURCE    | RESULT

        Integer     | Integer   | true
        Integer     | Float     | false
        Float       | Integer   | true
        Float       | Float     | true

        String      | String    | true
        String      | Integer   | false
        Integer     | String    | false

        Iterable    | Bag       | true
        Iterable    | List      | true
        Iterable    | Set       | true
    }

}
