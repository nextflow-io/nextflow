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

package nextflow.container

import groovy.transform.EqualsAndHashCode
import groovy.transform.ToString
import groovy.transform.TupleConstructor
/**
 * Hold container definition parsed in a shell script
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString
@EqualsAndHashCode
@TupleConstructor
class ContainerScriptTokens {

    static private final BASH_VAR = ~/([a-zA-Z_][0-9a-zA-Z_]*)=(.*)/

    String image

    int index

    Map<String,String> variables

    List<String> lines

    /**
     * Parse a shell script retuning a triple which contains the docker image,
     * the line on which the docker image has been defined and any environment variable
     * eventually defined before the image definition
     *
     * @param str A shell script to parse
     * @return
     */
    static ContainerScriptTokens parse(String str) {
        def norm = str.stripIndent().trim()
        if( !norm )
            throw new IllegalArgumentException("Missing `container` name in script definition")

        def index=-1
        def vars = [:]
        def image = null
        def lines = norm.readLines()
        for( int i=0; i<lines.size(); i++ ) {
            def entry = lines[i].trim()
            // check if it's a comment
            if( entry.startsWith('#') ) continue

            // check if it's a var definition
            def matcher = entry =~ BASH_VAR
            if( matcher.matches() ) {
                vars.put( matcher.group(1), removeQuotes(matcher.group(2)) )
                continue
            }

            if( entry.size() == 0 )
                continue

            // otherwise must be the container name
            image = entry.tokenize().get(0)
            index = i
            break
        }

        if( !image )
            throw new IllegalArgumentException("Missing `container` name in script definition")

        // returns the first string token on the first line
        return new ContainerScriptTokens(image, index, vars, lines)
    }

    static private String removeQuotes( String str ) {
        if( !str )
            return str

        if( str.startsWith("'") && str.endsWith("'") )
            return str.substring(1, str.size()-1)

        if( str.startsWith('"') && str.endsWith('"'))
            return str.substring(1, str.size()-1)

        return str
    }

}
