/*
 * Copyright 2013-2018, Centre for Genomic Regulation (CRG)
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

package nextflow.io

import groovy.transform.CompileStatic

/**
 * Implements a  {@link InputStream} object backed on a {@link DataInput} instance
 *
 * @author Paolo Di Tommaso <paolo.ditommaso@gmail.com>
 */
@CompileStatic
class DataInputStreamAdapter extends InputStream {

    final private DataInput target

    DataInputStreamAdapter( DataInput input ) { this.target = input }

    @Override
    int read() throws IOException {
        try {
            return target.readUnsignedByte()
        }
        catch( EOFException | IndexOutOfBoundsException e ) {
            return -1
        }
    }

    long skip(long n) {
        target.skipBytes(n as int)
    }
}
