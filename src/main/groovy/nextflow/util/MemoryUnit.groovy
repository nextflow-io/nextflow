/*
 * Copyright (c) 2013-2014, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2014, Paolo Di Tommaso and the respective authors.
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

import java.text.DecimalFormat

import groovy.transform.EqualsAndHashCode

/**
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@EqualsAndHashCode(includes = 'size')
class MemoryUnit implements Comparable<MemoryUnit> {

    static private final FORMAT = /([0-9\.]+)\s*(\S)?B?/

    static private final List UNITS = [ "B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB" ]

    long size

    MemoryUnit( long value ) {
        this.size = value
    }

    MemoryUnit( String str ) {

        def matcher = (str =~~ FORMAT)
        if( !matcher.matches() ) {
            throw new IllegalArgumentException("Not a valid FileSize value: '$str'")
        }

        final value = matcher[0][1]?.toString()
        final unit = matcher[0][2]?.toString()?.toUpperCase()

        if ( !unit || unit == "B" ) {
            size = Long.parseLong(value)
        }
        else {
            int p = UNITS.indexOf(unit)
            if ( p == -1 ) {
                // try adding a 'B' specified
                unit += 'B'
                p = UNITS.indexOf(unit)
                if( p == -1 ) {
                    throw new IllegalArgumentException("Not a valid file size unit: ${str}")
                }
            }

            Math.with {
                size = round( Double.parseDouble(value) * pow(1024, p) )
            }
        }

    }

    long toBytes() {
        size
    }


    def String toString() {
        if(size <= 0) {
            return "0"
        }

        int digitGroups = (int) (Math.log10(size) / Math.log10(1024))
        new DecimalFormat("#,##0.#").format(size/Math.pow(1024, digitGroups)) + " " + UNITS[digitGroups]
    }

    @Override
    int compareTo(MemoryUnit that) {
        return this.size <=> that.size
    }
}
