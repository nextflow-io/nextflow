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
/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
class CmdLineHelper {

    def List<String> args

    CmdLineHelper( String cmdLineToBeParsed ) {
        args = splitter(cmdLineToBeParsed ?: '')
    }

    def boolean contains(String argument) {
        return args.indexOf(argument) != -1
    }

    def getArg( String argument ) {
        def pos = args.indexOf(argument)
        if( pos == -1 ) return null

        def result = []
        for( int i=pos+1; i<args.size(); i++ ) {
            if( args[i].startsWith('-') ) {
                break
            }
            result.add(args[i])
        }

        if( result.size()==0 ) {
            return true
        }
        else if( result.size()==1 ) {
            return result[0]
        }
        else {
            return result
        }
    }

    def asList( String argument, String splitter=',' ) {
        def val = getArg(argument)
        if( !val ) return val

        if( val instanceof Boolean ) {
            return []
        }

        if( val instanceof String ) {
            val = [val]
        }

        for( int i=0; i<val.size(); i++ ) {
            val[i] = val[i] ?. split(splitter)
        }

        return val.flatten()
    }


    /**
     * Given a string the splitter method separate it by blank returning a list of string.
     * Tokens wrapped by a single quote or double quotes are considered as a contiguous string
     * and is added as a single element in the returned list.
     * <p>
     * For example the string: {@code "alpha beta 'delta gamma'"} will return the following result
     * {@code ["alpha", "beta", "delta gamma"]}
     *
     * @param cmdline The string to be splitted in single elements
     * @return A list of string on which each entry represent a command line argument, or an
     * empty list if the {@code cmdline} parameter is empty
     */
    static List<String> splitter( String cmdline ) {

        List<String> result = []

        if( cmdline ) {
            QuoteStringTokenizer tokenizer = new QuoteStringTokenizer(cmdline);
            while( tokenizer.hasNext() ) {
                result.add(tokenizer.next())
            }
        }

        result
    }


    static String toLine( String... args ) {
        toLine( args as List)
    }

    static String toLine( List<String> args ) {
        def result = new ArrayList(args.size())
        args.each { result << (it.contains(' ') ? "'${it}'" : it) }
        return result.join(' ')
    }

}
