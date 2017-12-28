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
