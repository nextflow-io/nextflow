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
import groovy.transform.InheritConstructors
import groovy.util.logging.Slf4j
/**
 * Splits a generic byte array in chunks having the specified length
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@Slf4j
@CompileStatic
@InheritConstructors
class BytesSplitter extends AbstractBinarySplitter {

    @Override
    def process( InputStream targetObject) {
        assert targetObject != null

        def result = null
        int c=0
        long bytesCount=0
        byte[] buffer = new byte[counter.size]
        int item

        try {

            while( (item=targetObject.read()) != -1 ) {
                buffer[c++] = (byte)item

                if ( counter.isChunckComplete() ) {
                    result = invokeEachClosure(closure, buffer)
                    buffer = new byte[counter.size]
                    c = 0
                    counter.reset()
                }

                // -- check the limit of allowed rows has been reached
                if( limit>0 && ++bytesCount == limit )
                    break
            }

        }
        finally {
            targetObject.closeQuietly()
        }


        /*
         * if there's something remaining in the buffer it's supposed
         * to be the last entry
         */

        if ( c ) {
            if( c != counter.size ) {
                def copy = new byte[c]
                System.arraycopy(buffer,0,copy,0,c)
                buffer = copy
            }

            result = invokeEachClosure(closure, buffer)
        }

        return result
    }

    @Override
    protected CollectorStrategy createCollector() {
        return null
    }
}
