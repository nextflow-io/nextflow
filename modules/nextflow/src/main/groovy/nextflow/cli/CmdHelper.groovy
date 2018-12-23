/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.cli

/**
 * Cmd CLI helper methods
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdHelper {

    static private String operators = '<>!='

    /**
     * Expand single `=` to `==` in filter string
     *
     * @param filter A filter string e.g. `cpus = 4`
     * @return A filter in which `=` char is expanded to `==` operator
     */
    static String fixEqualsOp( String filter ) {
        if( !filter ) return filter

        def result = new StringBuilder()
        int i=0
        int len=filter.length()

        while( i<len ) {
            def ch = filter[i++]
            result << ch
            if( i<len-1 && filter[i]=='=' && filter[i+1]!='=' && !operators.contains(ch)) {
                result.append('=')
            }
        }
        return result.toString()
    }

}
