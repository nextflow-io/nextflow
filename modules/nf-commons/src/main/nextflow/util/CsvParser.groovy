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
import groovy.transform.PackageScope
import groovy.transform.ToString
import org.apache.commons.lang.StringUtils
/**
 * Simple comma-separated-values parser
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@ToString(includeFields = true, includeNames = true)
@CompileStatic
class CsvParser {

    final static private char COMMA = ',' as char

    private char quote

    private char separator = COMMA

    private String empty

    private boolean strip

    CsvParser setQuote(char ch) {
        this.quote = ch
        return this
    }

    CsvParser setQuote(String ch) {
        this.quote = firstChar(ch)
        return this
    }

    CsvParser setSeparator(char ch) {
        this.separator = ch ?: COMMA
        return this
    }

    CsvParser setSeparator(String ch) {
        this.separator = firstChar(ch) ?: COMMA
        return this
    }

    CsvParser setStrip(boolean value) {
        this.strip = value
        return this
    }

    static private char firstChar(String str) {
        if( !str )
            return 0 as char
        if( str.size()>1 )
            throw new IllegalArgumentException("Not a valid CVS character: $str")
        str.charAt(0)
    }

    List<String> parse( String line ) {
        def result = []
        while( line != null ) {
            if( !line ) {
                result.add(empty)
                break
            }
            else if( quote && line.charAt(0)==quote ) {
                line = readQuotedValue(line, result)
            }
            else {
                line = readSimpleValue(line, result)
            }
        }
        return result
    }

    private String readSimpleValue(String line, List<String> result) {
        def p = line.indexOf( (int)separator )
        if( p == -1 ) {
            result.add(stripBlanks(line))
            return null
        }
        else {
            result.add(stripBlanks(line.substring(0,p)) ?: empty)
            return line.substring(p+1)
        }
    }

    private String readQuotedValue(String line, List<String> result) {
        def value = line.substring(1)
        def p = value.indexOf( (int)quote )
        if( p == -1 )
            throw new IllegalStateException("Missing double-quote termination in CSV value -- offending line: $line")
        result.add( stripBlanks(value.substring(0,p)) )

        def next = p+1
        if( next<value.size() ) {
            if( value.charAt(next) != separator )
                throw new IllegalStateException("Invalid CSV value -- offending line: $line")
            return value.substring(next+1)
        }
        else {
            return null
        }
    }

    @PackageScope String stripBlanks(String str) {
        final value = strip ? StringUtils.strip(str) : str

        if( !quote )
            return value

        if( value.length()>1 && value.charAt(0)==quote && value.charAt(value.size()-1)==quote )
            return value.substring(1, value.length()-1)

        else
            return value
    }

}
