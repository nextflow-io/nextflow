/*
 * Copyright (c) 2013-2015, Centre for Genomic Regulation (CRG).
 * Copyright (c) 2013-2015, Paolo Di Tommaso and the respective authors.
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

    private static char[] CHAR_ARRAY = new char[0]

    private LinkedList<String> holder = new LinkedList<String>()

    private long count

    @Override
    void add(Object record) {
        if( record==null )
            return

        def str = record.toString()
        count += str.length()
        holder.add(str)
    }

    @Override
    boolean hasChunk() {
        return !holder.isEmpty()
    }

    @Override
    void next() {
        count = 0
        holder.clear()
    }

    @Override
    def getChunk() {
        toString()
    }

    String toString() {

        if( count > Integer.MAX_VALUE )
            throw new IllegalArgumentException("String cannot contain more than ${Integer.MAX_VALUE} characters -- Your string needs ${count} characters")

        int offset = 0
        char[] content = new char[count]

        for( int i=0; i<holder.size(); i++ ) {
            final str = holder.get(i)
            final len = str.length()
            str.getChars(0,len, content, offset)
            offset += len
        }


        try {
            // use package private constructor to avoid to avoid to create a copy
            // of the buffer character array
            def constructor = String.class.getDeclaredConstructor(CHAR_ARRAY.class, boolean.class)
            constructor.setAccessible(true)
            constructor.newInstance(content, Boolean.TRUE)
        }
        catch( NoSuchMethodException e ) {
            // the above private constructor is not available is not available in all 1.7 runtime versions
            // fallback to the standard String constructor
            return new String(content)
        }
    }
}
