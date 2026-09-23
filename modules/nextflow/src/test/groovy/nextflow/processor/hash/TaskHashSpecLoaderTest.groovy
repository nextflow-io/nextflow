/*
 * Copyright 2013-2026, Seqera Labs
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
package nextflow.processor.hash

import spock.lang.Specification
import spock.lang.Unroll

class TaskHashSpecLoaderTest extends Specification {

    private static String spec(String keys, String function = 'murmur3_128') {
        return """{
          "id": "test/v1",
          "function": "${function}",
          "encoding": { "orderIndependentMaps": true, "cacheFunnelFirst": true },
          "keys": ${keys}
        }"""
    }

    def 'loads a well-formed spec'() {
        when:
        def result = TaskHashSpecLoader.load(spec('''[
            {"key":"SESSION_ID","contributor":"sessionId"},
            {"key":"CONDA","contributor":"condaEnv"}]'''), 'test')

        then:
        result.id == 'test/v1'
        result.keys() == [HashKey.SESSION_ID, HashKey.CONDA]
        result.encoding.canonicalForm() == 'orderIndependentMaps=true;cacheFunnelFirst=true'
    }

    @Unroll
    def 'rejects #reason rather than loading a spec nobody wrote'() {
        when:
        TaskHashSpecLoader.load(json, 'test')
        then:
        thrown(IllegalArgumentException)

        where:
        reason                  | json
        'an unknown contributor'  | spec('[{"key":"SESSION_ID","contributor":"nope"}]')
        'an unknown hash key'   | spec('[{"key":"NOT_A_KEY","contributor":"sessionId"}]')
        'a duplicated key'      | spec('''[{"key":"CONDA","contributor":"condaEnv"},
                                           {"key":"CONDA","contributor":"condaEnv"}]''')
        'an entry missing its contributor' | spec('[{"key":"CONDA"}]')
        'an unsupported function'| spec('[{"key":"CONDA","contributor":"condaEnv"}]', 'sha256')
        'empty keys'            | spec('[]')
        'a non-object document' | '[]'
    }

    def 'every shipped std spec loads and every contributor it names is registered'() {
        expect:
        StdSpecs.all().size() == 4
        and:
        StdSpecs.all().every { s -> s.bindings.every { b -> b.contributor.canonicalName() in ContributorRegistry.names() } }
    }
}
