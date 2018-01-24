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
