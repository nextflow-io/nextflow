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
