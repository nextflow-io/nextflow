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

package nextflow.splitter

import groovy.transform.CompileStatic

/**
 * A collector strategy that creates a string chunks
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CharSequenceCollector implements CollectorStrategy {

    private static final int INITIAL_SIZE = 10 * 1024 * 1024

    private StringBuilder holder = new StringBuilder(INITIAL_SIZE)

    private long count

    @Override
    void add(Object record) {
        if( record==null )
            return

        def str = record.toString()
        count += str.length()
        holder.append(str)
    }

    @Override
    boolean hasChunk() {
        return holder.length()>0
    }

    @Override
    def nextChunk() {
        try {
            return toString()
        }
        finally {
            count = 0
            holder.setLength(0)
        }
    }

    String toString() {

        if( count > Integer.MAX_VALUE )
            throw new IllegalArgumentException("String cannot contain more than ${Integer.MAX_VALUE} characters -- Your string needs ${count} characters")

        holder.toString()
    }
}
