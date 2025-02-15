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
package nextflow.plugin

import org.pf4j.PluginManager
import org.pf4j.PluginWrapper
import spock.lang.Specification
import spock.lang.Unroll

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class BasePluginTest extends Specification {

    @Unroll
    def 'should validate version match' () {
        given:
        def wrapper = Mock(PluginWrapper) {
            getPluginManager() >> Mock(PluginManager) {
                getVersionManager() >> new CustomVersionManager()
            }
        }

        def plugin = new BasePlugin(wrapper) { }

        expect:
        plugin.verMatches(REQUIRES, CURRENT) == EXPECT

        where:
        CURRENT         | REQUIRES                      | EXPECT
        '20.01.0'       | '21.01.0'                     | false
        '21.01.0'       | '21.01.0'                     | true
        '21.01.0'       | '>=21.01.0'                   | true
        '21.01.0-edge'  | '>=21.01.0-edge'              | true
        '21.01.0-edge'  | '*'                           | true
        and:
        '21.01.0-edge'  | 'nextflow@>=21.01.0-edge'     | true
        '20.01.0-edge'  | 'nextflow@>=21.01.0-edge'     | false
        '21.01.0-SNAPSHOT'  | '>=21.01.0-edge'          | true
    }
}
