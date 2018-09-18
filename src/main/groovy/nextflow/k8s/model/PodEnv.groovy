/*
 * Copyright (c) 2013-2018, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2018, Paolo Di Tommaso and the respective authors.
 *
 *   This file is part of 'Nextflow'.
 *
 *   Nextflow is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   Nextflow is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with Nextflow.  If not, see <http://www.gnu.org/licenses/>.
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
