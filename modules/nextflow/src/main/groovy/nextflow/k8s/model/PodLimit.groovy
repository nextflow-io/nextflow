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
 * Model a Pod limits
 * @author Kevin Sayers <sayerskt@gmail.com>
 */

@CompileStatic
@ToString(includeNames = true)
@EqualsAndHashCode(includeFields = true)
class PodLimit {

    private Map spec = [:]

    PodLimit(limit) {
        if( limit instanceof CharSequence )
            createWithString(limit.toString())

        else if( limit instanceof Map )
            createWithMap(limit)

        else if( limit != null )
            throw new IllegalArgumentException("Not a valid pod limit value: $limit [${limit.getClass().getName()}]")
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
    private createWithString( String limit ) {
        if(!limit) return
        def entries = limit.tokenize(',')
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
        "PodLimit[ ${spec?.toString()} ]"
    }

}


