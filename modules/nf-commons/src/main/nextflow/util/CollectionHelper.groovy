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

package nextflow.util

import groovy.transform.CompileStatic

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CollectionHelper {


    static List flatten( Collection collection ) {
        def result = []
        flatten(collection, { result << it })
        return result
    }

    static void flatten( Collection list, Closure row ) {

        def maxLen = 1
        list.each { if( it instanceof Collection ) { maxLen = Math.max(maxLen,((Collection)it).size()) } }

        for( int i=0; i<maxLen; i++ ) {
            def set = new ArrayList(list.size())

            for( int j=0; j<list.size(); j++ ) {
                def val = list[j]
                if( val instanceof Collection )
                    set[j] = i < ((Collection)val).size() ? ((Collection)val)[i] : null
                else
                    set[j] = val
            }

            row.call(set)
        }


    }

}
