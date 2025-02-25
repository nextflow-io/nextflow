/*
 * Copyright 2013-2024, Seqera Labs
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

package nextflow.secret

import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class SecretHelperTest extends Specification {

    @Unroll
    def 'should validate name #NAME'() {
        when:
        SecretsHelper.checkName(NAME)
        then:
        noExceptionThrown()

        where:
        _ | NAME
        _ | 'ab'
        _ | 'AB'
        _ | 'A123'
        _ | 'A_123'
        _ | '_A123'
        _ | 'ABC_123'
        _ | 'ABC_123_'
    }

    @Unroll
    def 'should not validate name #NAME'() {
        when:
        SecretsHelper.checkName(NAME)
        then:
        thrown(IllegalArgumentException)

        where:
        _ | NAME
        _ | 'a'
        _ | '0'
        _ | '0a'
        _ | 'A__B'
        _ | 'AB__'
        _ | '__AB'

    }


}
