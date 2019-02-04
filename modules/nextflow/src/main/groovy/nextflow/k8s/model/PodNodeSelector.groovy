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
 * Model a Pod nodeSelector spec
 *
 *  https://kubernetes.io/docs/concepts/configuration/assign-pod-node/#nodeselector
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode(includeFields = true)
class PodNodeSelector {

    private Map spec = [:]

    PodNodeSelector(selector) {
        if( selector instanceof CharSequence )
            createWithString(selector.toString())

        else if( selector instanceof Map )
            createWithMap(selector)

        else if( selector != null )
            throw new IllegalArgumentException("Not a valid pod nodeSelector value: $selector [${selector.getClass().getName()}]")
    }

    private createWithMap(Map selection ) {
        if(selection) {
            for( Map.Entry entry : selection.entrySet() ) {
                spec.put(entry.key.toString(), entry.value?.toString())
            }
        }
    }

    /**
     * @param selector
     *      A string representing a comma separated list of pairs
     *      e.g. foo=1,bar=2
     *
     */
    private createWithString( String selector ) {
        if(!selector) return
        def entries = selector.tokenize(',')
        for( String item : entries ) {
            def pair = item.tokenize('=')
            spec.put( trim(pair[0]), trim(pair[1]) ?: 'true' )
        }
    }

    private String trim(String v) {
        v?.trim()
    }

    Map<String,String> toSpec() { spec }

    String toString() {
        "PodNodeSelector[ ${spec?.toString()} ]"
    }
}
