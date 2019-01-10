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

package nextflow.splitter

import groovy.transform.CompileStatic

/**
 * A collector strategy that creates a string chunks
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class CharSequenceCollector implements CollectorStrategy, HeaderCollector {

    private static final int INITIAL_SIZE = 10 * 1024 * 1024

    private StringBuilder holder = new StringBuilder(INITIAL_SIZE)

    private long count

    private String header

    void setHeader(String value) {
        this.header = value
    }

    @Override
    void add(Object record) {
        if( record==null )
            return

        if( count==0 && header ) {
            holder.append(header)
        }
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
