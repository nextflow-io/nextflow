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

package nextflow.k8s.model

import groovy.transform.CompileStatic
import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString

/**
 * Model a K8s pod environment variable definition
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode(includeFields = true)
class PodEnv {

    private Map spec

    private PodEnv(Map spec) {
        this.spec = spec
    }

    static PodEnv value(String env, String value) {
        new PodEnv([name:env, value:value])
    }

    static PodEnv config(String env, String config) {
        final tokens = config.tokenize('/')
        if( tokens.size() > 2 )
            throw new IllegalArgumentException("K8s invalid pod env file: $config -- Secret must be specified as <config-name>/<config-key>")

        final name = tokens[0]
        final key = tokens[1]

        assert env, 'Missing pod env variable name'
        assert name, 'Missing pod env config name'

        final ref = [ name: name, key: (key ?: env) ]
        new PodEnv([ name: env, valueFrom: [configMapKeyRef: ref]])
    }

    static PodEnv secret(String env, String secret) {

        final tokens = secret.tokenize('/')
        if( tokens.size() > 2 )
            throw new IllegalArgumentException("K8s invalid pod env secret: $secret -- Secret must be specified as <secret-name>/<secret-key>")

        final name = tokens[0]
        final key = tokens[1]

        final ref = [ name: name, key: (key ?: env) ]
        new PodEnv([ name: env, valueFrom: [secretKeyRef: ref]])
    }


    Map toSpec() { spec }

    String toString() {
        "PodEnv[ ${spec?.toString()} ]"
    }
}
