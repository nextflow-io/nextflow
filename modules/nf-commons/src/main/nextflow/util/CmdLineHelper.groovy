/*
 * Copyright 2013-2024, Seqera Labs
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

import java.util.regex.Pattern

/**
 * Implement command line parsing helpers
 * 
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CmdLineHelper {

    static private Pattern CLI_OPT = ~/--([a-zA-Z_-]+)(?:\W.*)?$|-([a-zA-Z])(?:\W.*)?$/

    private List<String> args

    CmdLineHelper( String cmdLineToBeParsed ) {
        args = splitter(cmdLineToBeParsed ?: '')
    }

    private boolean contains(String argument) {
        return args.indexOf(argument) != -1
    }

    private getArg( String argument ) {
        int pos = args.indexOf(argument)
        if( pos == -1 ) return null

        List<String> result = []
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

    /**
     * Parse command line and returns the options and their values as a map object.
     *
     * @param cmdline
     *      The command line as single string
     * @return
     *      A map object holding the option key-value(s) associations
     */
    static CmdLineOptionMap parseGnuArgs(String cmdline) {
        final BLANK = ' ' as char
        final result = new CmdLineOptionMap()

        if( !cmdline )
            return result

        final tokenizer = new QuoteStringTokenizer(cmdline, BLANK);
        String opt = null
        String last = null
        while( tokenizer.hasNext() ) {
            final String token = tokenizer.next()
            if( !token || token=='--')
                continue
            final matcher = CLI_OPT.matcher(token)
            if(  matcher.matches() ) {
                if( opt ) {
                    result.addOption(opt,'true')
                }
                opt = matcher.group(1) ?: matcher.group(2)
            }
            else {
                if( !opt ) {
                    if( !last ) continue
                    result.addOption(last, token)
                }
                else {
                    result.addOption(opt, token)
                    last = opt
                    opt = null
                }
            }
        }

        if( opt )
            result.addOption(opt, 'true')

        return result
    }
}
